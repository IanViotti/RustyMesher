use crate::config;
use ndarray::Array2;
use std::f64::consts::PI;
use crate::mesher_utils::Point;
use crate::config::AirfoilType;

/// Main entry point for setting up the physical boundaries of the computational domain.
/// This applies the Dirichlet boundary conditions for the grid generation:
/// 1. Inner boundary (j = 0): The airfoil surface.
/// 2. Outer boundary (j = j_max): The far-field circular boundary.
pub fn insert_geoemtry(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Inserting geometry...");

    // Routes to the correct inner boundary generation function based on the selected enum.
    // The inner boundary is located at the lowest normal index (eta / j = 0).
    match config.airfoil_type {
        AirfoilType::Biconvex => insert_biconvex_geom(config, grid),
        AirfoilType::NACA00XX => insert_naca00xx_geom(config, grid),
    }

    // Insert the outer boundary geometry into the grid (eta / j = j_max)
    insert_outter_boundary(config, grid);
}

/// Generates the coordinates for a symmetrical biconvex airfoil.
/// The profile is mathematically defined as a parabola: y = 2 * t * x * (1 - x)
fn insert_biconvex_geom(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Generating Biconvex airfoil geometry...");

    let n_longitudinal = config.longitudinal_points;
    let t = config.t;
    
    // The midpoint represents the Leading Edge (LE) in the O-grid topology.
    // i = 0 is the Trailing Edge (TE) lower surface, and i = n_longitudinal-1 is the TE upper surface.
    let mid = (n_longitudinal - 1) / 2;
    
    // Call the helper function to get the clustered x-coordinate distribution [0.0 to 1.0]
    let x_dist = generate_x_distribution(config);

    for i in 0..n_longitudinal {
        // Map the O-grid longitudinal index 'i' to the structural coordinate 'x'.
        // It wraps around the airfoil: from TE (x=1) -> LE (x=0) -> TE (x=1).
        let x = if i <= mid { x_dist[mid - i] } else { x_dist[i - mid] };
        let x_safe = x.max(0.0);

        // Biconvex thickness equation
        let thickness = 2.0 * t * x_safe * (1.0 - x_safe);

        // Assign negative thickness for the lower surface (i <= mid) and positive for the upper (i > mid)
        let y = if i <= mid { -thickness } else { thickness };

        grid[[i, 0]] = Point { x, y };
    }
}

/// Generates the coordinates for a standard NACA 4-digit symmetrical airfoil (e.g., NACA 0012).
fn insert_naca00xx_geom(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Generating NACA 00XX");

    let n_longitudinal = config.longitudinal_points;
    let t = config.t;
    let mid = (n_longitudinal - 1) / 2;
    
    // Vector with the mathematical point distribution (from 0.0 to 1.0)
    let x_dist = generate_x_distribution(config);

    for i in 0..n_longitudinal {
        // x_math dictates the position along the chord.
        // 1.0 is the trailing edge (tail), 0.0 is the leading edge (nose).
        let x_math = if i <= mid {
            x_dist[mid - i] // Descending from TE (1.0) to LE (0.0) along the lower surface
        } else {
            x_dist[i - mid] // Ascending from LE (0.0) back to TE (1.0) along the upper surface
        };

        let x_safe = x_math.max(0.0);

        // Standard NACA 4-digit thickness equation.
        // The last coefficient (-0.1036) ensures a closed/sharp trailing edge at x = 1.0.
        // (If -0.1015 was used, it would leave a finite trailing edge thickness).
        let thickness = 5.0 * t * (
            0.2969 * x_safe.sqrt() - 
            0.1260 * x_safe - 
            0.3516 * x_safe.powi(2) + 
            0.2843 * x_safe.powi(3) - 
            0.1036 * x_safe.powi(4)
        );

        // === THE ONLY CHANGE TO POINT TO THE LEFT ===
        // When x_math is 1.0 (tail), physical X will be 1.0 (or 0.0 depending on origin).
        // Kept as x_safe to match the coordinate system orientation.
        let x_coord = x_safe;

        // Strictly maintains the sweep direction required by O-grid topology:
        // i <= mid: Goes along the lower surface (pressure side / intrados)
        // i > mid: Returns along the upper surface (suction side / extrados)
        let y_coord = if i <= mid {
            -thickness 
        } else {
            thickness  
        };

        grid[[i, 0]] = Point { x: x_coord, y: y_coord };
    }
}

/// Generates the normalized X-coordinate distribution [0.0 to 1.0] with algebraic stretching.
/// Clustering points near the Leading Edge (x=0) and Trailing Edge (x=1) is crucial in CFD 
/// to accurately capture high flow gradients (e.g., stagnation points, shocks, wakes).
fn generate_x_distribution(config: &config::Config) -> Vec<f64> {
    let n_longitudinal = config.longitudinal_points;
    let r = config.stretching_factor; 

    // We only need to calculate half of the chord (0.0 to 1.0), then mirror it later.
    let mid = (n_longitudinal - 1) / 2;
    let num_pts = mid + 1; 
    let mut x_dist = vec![0.0; num_pts];

    if r <= 1.0 {
        // Uniform grid distribution (No clustering)
        for k in 0..num_pts {
            x_dist[k] = (k as f64) / ((num_pts - 1) as f64);
        }
    } else {
        // Stretching using Geometric Progression (Clustered at edges)
        let n_half = (num_pts - 1) / 2; 
        
        // Calculate the initial step size (a1) ensuring the sum of the geometric series reaches 0.5 at mid-chord
        let a1 = 0.5 * (r - 1.0) / (r.powi(n_half as i32) - 1.0);
        
        x_dist[0] = 0.0;
        let mut current_dx = a1;
        
        // First quarter: from LE (x=0.0) to mid-chord (x=0.5)
        for k in 1..=n_half {
            x_dist[k] = x_dist[k - 1] + current_dx;
            current_dx *= r; 
        }
        
        // Second quarter: Mirror the distribution from mid-chord (x=0.5) to TE (x=1.0)
        for k in (n_half + 1)..num_pts {
            let mirror_idx = (num_pts - 1) - k;
            x_dist[k] = 1.0 - x_dist[mirror_idx];
        }
    }

    x_dist
}

/// Defines the far-field (outer boundary) geometry.
/// For an O-grid, this is typically a large circle centered around the airfoil.
fn insert_outter_boundary(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Inserting external boundary...");

    let n_longitudinal = config.longitudinal_points;
    
    // Maximum radius of the computational domain (distance to far-field)
    let r_max = config.r_max; 

    // Obtain the normal dimension size (eta) directly from the grid matrix shape
    let (_, n_normal) = grid.dim();
    let j_max = n_normal - 1; // Index of the outermost boundary layer

    // Center of the external circle (typically placed at mid-chord for symmetrical grids)
    let x_center = 0.5;
    let y_center = 0.0;

    for i in 0..n_longitudinal {
        // Angular sweep from 0 to -2*PI.
        // It is CRITICAL that the direction of this sweep matches the inner boundary loop.
        // Sweeping in the negative direction ensures proper cell orientation (positive Jacobian).
        // i = 0 (Lower TE)       -> angle 0
        // i = mid (LE)           -> angle -PI
        // i = final (Upper TE)   -> angle -2*PI
        let theta = -2.0 * PI * (i as f64) / ((n_longitudinal - 1) as f64);

        // Parametric equation of a circle
        let x = x_center + r_max * theta.cos();
        let y = y_center + r_max * theta.sin();

        // Assign the computed far-field points to the outermost row of the grid (j_max)
        grid[[i, j_max]] = Point { x, y };
    }
}