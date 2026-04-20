use std::vec;
use crate::config::Config;
use ndarray::Array2;
use crate::mesher_utils::Point;
use crate::mesher_utils::periodic_thomas;

/// Entry point for generating the initial grid using a parabolic marching method.
/// This method generates the internal grid points algebraically/parabolically 
/// by marching from the inner boundary (airfoil) to the outer boundary.
/// It provides a high-quality, non-overlapping initial guess for the elliptic smoother.
pub fn generate_parabolic_mesh(config: &Config, grid: &mut Array2<Point>) {
    println!("Generating parabolic mesh...");
    
    space_march(config, grid);

}

/// Performs the space-marching algorithm (typically based on Nakamura's approach).
/// It solves a set of parabolic partial differential equations layer by layer 
/// (from j=1 to j_max), advancing outwards from the airfoil surface.
fn space_march(config: &Config, grid: &mut Array2<Point>) {
    let ni = config.longitudinal_points - 1;
    let nj = config.normal_points - 1;

    // Retrieve the normalized stretching distribution for the normal (eta) direction
    let s_bar = streching_factor(config);

    // March outwards layer by layer
    for j in 1..nj{

        // Allocate arrays for the tridiagonal system coefficients (a, b, c) 
        // and right-hand side vectors (d_x, d_y)
        //let mut ref_grid_jp1: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];
        //let mut ref_grid_j: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];
        let mut a = vec![0.0; ni];
        let mut b = vec![0.0; ni];
        let mut c = vec![0.0; ni];
        let mut d_x = vec![0.0; ni];
        let mut d_y = vec![0.0; ni];

        // Grid spacing control in eta direction
        //let ds = (s_bar[j] - s_bar[j-1]) / (s_bar[nj] - s_bar[j-1]);

        // Calculate the local reference grid points for the current and next layer (j and j+1).
        // This guides the marching scheme, ensuring lines blend smoothly towards the outer boundary.
        let (ref_grid_j, ref_grid_jp1) = local_reference_grid(config, grid, j, &s_bar);

        // Loop over the longitudinal direction (wrap-around the airfoil)
        for i in 0..ni {

            // Define adjacent indices (i-1, i+1) handling the periodic boundary condition
            // of the O-grid topology where the wake cut connects the start and end of the loop.
            let ip1: usize;
            let im1: usize;

            if i == 0 {
                ip1 = 1;
                im1 = ni - 1;
            } else if i == ni - 1{
                ip1 = 0;
                im1 = i - 1;
            }
            else{
                ip1 = i + 1;
                im1 = i - 1;
            }

            // Central differences for derivatives with respect to xi (longitudinal)
            let x_xi = 1.0 / 2.0 * (ref_grid_jp1[ip1].x - ref_grid_jp1[im1].x);
            let y_xi = 1.0 / 2.0 * (ref_grid_jp1[ip1].y - ref_grid_jp1[im1].y);

            // Forward differences for derivatives with respect to eta (normal)
            let x_eta = ref_grid_jp1[i].x - ref_grid_j[i].x;
            let y_eta = ref_grid_jp1[i].y - ref_grid_j[i].y;

            // Compute the metric tensor coefficients (A, B, C) derived from the mapping Jacobian.
            // A controls eta diffusion, B controls cross-derivatives, C controls xi diffusion.
            let A = x_eta.powi(2) + y_eta.powi(2);
            let B = x_xi * x_eta + y_xi * y_eta;
            let C = x_xi.powi(2) + y_xi.powi(2);

            // Populate the coefficients for the implicit tridiagonal system
            a[i] = 2.0 * A;
            b[i] = -4.0 * (A + C);
            c[i] = 2.0 * A;
            
            // Right-hand side of the tridiagonal system (forcing terms)
            d_x[i] = B * (ref_grid_jp1[ip1].x - ref_grid_jp1[im1].x - grid[[ip1, j-1]].x + grid[[im1, j-1]].x) -
                        2.0 * C * (ref_grid_jp1[i].x + grid[[i, j-1]].x); 
            d_y[i] = B * (ref_grid_jp1[ip1].y - ref_grid_jp1[im1].y - grid[[ip1, j-1]].y + grid[[im1, j-1]].y) -
                        2.0 * C * (ref_grid_jp1[i].y + grid[[i, j-1]].y); 
        }

        // Solve the periodic tridiagonal system to obtain the physical coordinates x and y 
        // for the current marching layer j.
        let x = periodic_thomas(&a, &b, &c, &d_x);
        let y = periodic_thomas(&a, &b, &c, &d_y);

        // Update the main grid with the newly solved coordinates
        for i in 0..ni{
            grid[[i, j]].x = x[i]; 
            grid[[i, j]].y = y[i];

            //grid[[i, j]] = ref_grid_jp1[i]; 

        }
        // Overlap points at the wake cut to strictly enforce periodicity (closing the O-grid)
        grid[[ni, j]] = grid[[0, j]];

    }
}


/// Calculates the local reference grid for layers j and j+1.
/// The reference grid smoothly transitions the algebraic grid generation from 
/// being orthogonal to the airfoil near the inner boundary to interpolating 
/// radially towards the outer boundary, preventing grid lines from crossing.
/// Returns a tuple containing the point vectors for layers j and j+1.
fn local_reference_grid(config: &Config, grid: &mut Array2<Point>, j: usize, s_bar: &Vec<f64>) -> (Vec<Point>, Vec<Point>) {

    // Outer boundary index
    let nj = config.normal_points - 1;
    let ni = config.longitudinal_points - 1;

    // Allocate vectors for the reference layers 
    let mut ref_grid_j: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];
    let mut ref_grid_jp1: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];

    for i in 0..ni{
    
        // Periodic boundary indices
        let ip1: usize ;
        let im1: usize;

        if i == 0 {
            im1 = ni - 1;
            ip1 = i + 1;
        }
        else if i == ni - 1 {
            im1 = i - 1;
            ip1 = 0;
        }
        else {
            im1 = i - 1;
            ip1 = i + 1;
        }

        // === Local Reference Grid @ j ===

        // xi derivatives (central difference)
        let x_xi = 1.0 / 2.0 * (grid[[ip1, j-1]].x - grid[[im1, j-1]].x);
        let y_xi = 1.0 / 2.0 * (grid[[ip1, j-1]].y - grid[[im1, j-1]].y);

        // Vector pointing from the previous layer to the outer boundary
        let dx = grid[[i, nj]].x - grid[[i, j-1]].x;
        let dy = grid[[i, nj]].y - grid[[i, j-1]].y;
        
        // Local normalized spacing step
        let ds = (s_bar[j] - s_bar[j-1]) / (s_bar[nj] - s_bar[j-1]);
        
        // Arc-length constraint to enforce orthogonality at the wall
        let S_eta = calc_S_eta(dx, dy, ds);

        // eta derivatives representing a vector strictly normal to the xi-line
        let x_eta = - (y_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();
        let y_eta = (x_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();

        // Orthogonally projected coordinates
        let x_o = grid[[i, j-1]].x + x_eta;
        let y_o = grid[[i, j-1]].y + y_eta;

        // Straight-line interpolated coordinates to the outer boundary
        let x_int = grid[[i, j-1]].x + ds * (grid[[i, nj]].x - grid[[i, j-1]].x);
        let y_int = grid[[i, j-1]].y + ds * (grid[[i, nj]].y - grid[[i, j-1]].y);

        // Blending factor (eps): 0.0 near the airfoil (fully orthogonal), 
        // approaching 1.0 near the outer boundary (fully interpolated)
        let eps = (j as f64 - 1.0) / (nj as f64 - 1.0);

        // Compute the blended reference grid point
        ref_grid_j[i].x = eps * x_int + (1.0 - eps) * x_o;
        ref_grid_j[i].y = eps * y_int + (1.0 - eps) * y_o;

    }


    // === Local Reference Grid @ j + 1 ===

    for i in 0..ni{

        // Periodic boundary indices
        let ip1: usize ;
        let im1: usize;

        if i == 0 {
            im1 = ni - 1;
            ip1 = i + 1;
        }
        else if i == ni - 1 {
            im1 = i - 1;
            ip1 = 0;
        }
        else {
            im1 = i - 1;
            ip1 = i + 1;
        }

        // xi derivatives using the previously computed reference layer j
        let x_xi = 1.0 / 2.0 * (ref_grid_j[ip1].x - ref_grid_j[im1].x);
        let y_xi = 1.0 / 2.0 * (ref_grid_j[ip1].y - ref_grid_j[im1].y);

        // Vector pointing from reference layer j to the outer boundary
        let dx = grid[[i, nj]].x - ref_grid_j[i].x;
        let dy = grid[[i, nj]].y - ref_grid_j[i].y;
        
        let ds = (s_bar[j+1] - s_bar[j]) / (s_bar[nj] - s_bar[j]);
        let S_eta = calc_S_eta(dx, dy, ds);

        // Normal vector projection
        let x_eta = - (y_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();
        let y_eta = (x_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();

        // Orthogonally projected coordinates
        let x_o = ref_grid_j[i].x + x_eta;
        let y_o = ref_grid_j[i].y + y_eta;

        // Straight-line interpolated coordinates
        let x_int = ref_grid_j[i].x + ds * (grid[[i, nj]].x - ref_grid_j[i].x);
        let y_int = ref_grid_j[i].y + ds * (grid[[i, nj]].y - ref_grid_j[i].y);

        // Blending factor (eps) for layer j+1
        let eps = (j as f64 - 0.0) / (nj as f64 - 1.0);

        // Compute the blended reference grid point for j+1
        ref_grid_jp1[i].x = eps * x_int + (1.0 - eps) * x_o;
        ref_grid_jp1[i].y = eps * y_int + (1.0 - eps) * y_o;
    }

    (ref_grid_j, ref_grid_jp1)

}

/// Generates a geometric stretching factor vector for the normal direction (eta).
/// This controls how tightly packed the grid lines are near the airfoil surface.
fn streching_factor(config: &Config) -> Vec<f64> {
    let nj = config.normal_points - 1;
    let q = config.q; // Geometric progression ratio

    let mut s_bar = vec![0.0; nj + 1];
    
    // Index 0 represents the airfoil surface (s = 0.0)
    let mut current_ds = 1.0; // Initial unnormalized step size

    // Fill the distribution from index 1 to nj
    for j in 0..nj {
        s_bar[j+1] = s_bar[j] + current_ds;
        current_ds *= q; // The next step will be scaled by factor q
    }

    s_bar
}

/// Calculates the constraint scalar S_eta, which enforces the spacing and 
/// orthogonality distance for the reference grid algebraic projection.
fn calc_S_eta(dx: f64, dy: f64, ds: f64) -> f64 {

    // Distance vector magnitude
    let R = (dx.powi(2) + dy.powi(2)).sqrt();

    // The desired normal step size
    let S_eta = ds * R;

    S_eta
}