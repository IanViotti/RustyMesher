use crate::config;
use ndarray::Array2;
use std::fs::File;
use std::io::{self, Write};

/// Basic structure representing a 2D Cartesian coordinate in the physical domain.
#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

/// Initializes the output directory for the current mesh job.
/// Ensures that the required folder structure exists before saving files.
pub fn init_mesh_directory(config: &config::Config) {
    let meshname = &config.meshname;
    std::fs::create_dir_all(&format!("job_files/{}", meshname)).unwrap();
}

/// Allocates the 2D array representing the computational grid (xi, eta).
/// Initializes all points at the origin (0.0, 0.0) before the solver maps them to physical space.
pub fn create_grid(longitudinal_points: usize, normal_points: usize) -> Array2<Point> {
    // Initializes the matrix with dimensions (xi, eta), filling it entirely with Point {x: 0.0, y: 0.0}
    Array2::from_elem((longitudinal_points, normal_points), Point { x: 0.0, y: 0.0 })
}

/// Exports the grid coordinates to a basic CSV format.
/// Useful for quick data inspection, debugging, or custom Python plotting scripts.
pub fn save_grid_to_csv(grid: &Array2<Point>, filepath: &str) -> io::Result<()> {
    // Creates (or overwrites) the file at the specified path
    let mut file = File::create(filepath)?;

    // Writes the header manually
    writeln!(file, "i,j,x,y")?;

    // Iterates over all elements and their indices, writing line by line
    for ((i, j), point) in grid.indexed_iter() {
        writeln!(file, "{},{},{},{}", i, j, point.x, point.y)?;
    }

    Ok(())
}

/// Thomas algorithm (Tridiagonal Matrix Algorithm - TDMA)
///
/// Solves a tridiagonal linear system of equations efficiently in O(N) time:
///
///     a_j x_{j-1} + b_j x_j + c_j x_{j+1} = d_j
///
/// In CFD mesh generation, this is heavily used in implicit schemes like 
/// ADI (Alternating Direction Implicit) or SLOR (Successive Line Over-Relaxation)
/// to solve for coordinate corrections along a grid line (constant xi or eta).
pub fn thomas_algorithm(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let mut c_star = vec![0.0; n];
    let mut d_star = vec![0.0; n];
    let mut x = vec![0.0; n];

    // Step 1: Forward sweep - Implicitly building L and U
    // Modifies coefficients to eliminate the lower diagonal 'a'
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for i in 1..n {
        let m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        if i < n - 1 {
            c_star[i] = c[i] * m;
        }
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    // Step 2: Backward substitution
    // Solves for 'x' from the bottom up
    x[n - 1] = d_star[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    x
}

/// Solves a cyclic (periodic) tridiagonal system using the Sherman-Morrison formula.
/// 
/// This is necessary for O-grid topologies. In an O-grid, the wrap-around 
/// direction (longitudinal/xi) forms a closed loop, meaning the first node and the 
/// last node on a given constant-eta line are adjacent (or identical at the wake cut). 
/// Standard TDMA fails here; the Sherman-Morrison approach handles the cyclic connections.
pub fn periodic_thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    
    // Defines gamma (usually -b[0] to ensure numerical stability)
    let gamma = -b[0]; 
    
    // Creates the modified main diagonal (Matrix T)
    let mut bb = b.to_vec();
    bb[0] = b[0] - gamma;
    bb[n - 1] = b[n - 1] - a[0] * c[n - 1] / gamma;

    // Assembles vector u (only the endpoints contain values)
    //let mut u = vec![0.0; n];
    //u[0] = gamma;
    //u[n - 1] = a[0]; 
    // Assembles vector u (only the endpoints contain values)
    let mut u = vec![0.0; n];
    u[0] = gamma;
    u[n - 1] = c[n - 1]; // CORRECTION: It was a[0]

    // Vector v is modeled implicitly to save memory allocation:
    // v = [1.0, 0.0, ..., 0.0, c[n-1]/gamma]

    // Double Resolution:
    // 1. Solves Ty = d using the modified main diagonal
    let y: Vec<f64> = thomas_algorithm(a, &bb, c, d);
    
    // 2. Solves Tq = u
    let q = thomas_algorithm(a, &bb, c, &u);

    // 3. Final reconstruction of the solution combining the results
    //let v_n_minus_1 = c[n - 1] / gamma;
    // 3. Final reconstruction of the solution combining the results
    let v_n_minus_1 = a[0] / gamma; // CORRECTION: It was c[n - 1] / gamma
    
    // Dot products (v dot y) and (v dot q)
    let v_dot_y = y[0] + v_n_minus_1 * y[n - 1];
    let v_dot_q = q[0] + v_n_minus_1 * q[n - 1];

    let factor = v_dot_y / (1.0 + v_dot_q);

    let mut x = vec![0.0; n];
    for i in 0..n {
        x[i] = y[i] - factor * q[i];
    }

    x
}    

/// Exports the 2D grid to the VTK Legacy format (Structured Grid).
/// This is the standard output format to visualize the mesh topology, 
/// orthogonality, and quality using post-processing software like ParaView.
pub fn export_vtk_structured_grid(grid: &Array2<Point>, filename: &str) -> io::Result<()> {

    // Extracts the actual dimensions of the allocated mesh
    let shape = grid.shape();
    let ni = shape[0];
    let nj = shape[1];
    let num_points = ni * nj;

    // Creates the file (overwrites if it already exists)
    let mut file = File::create(filename)?;

    // --- 1. VTK Legacy Header ---
    // Standard VTK headers defining version, title, and data type
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "Malha Parabolica 2D O-Mesh")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_GRID")?;
    
    // VTK always expects 3 dimensions (i, j, k). For 2D, k=1.
    writeln!(file, "DIMENSIONS {} {} 1", ni, nj)?;
    writeln!(file, "POINTS {} float", num_points)?;

    // --- 2. Writing Coordinates ---
    // The order REQUIRES that 'i' be the innermost loop.
    // VTK requires points to be listed varying 'i' fastest, then 'j', then 'k'.
    for j in 0..nj {
        for i in 0..ni {
            let pt = &grid[[i, j]];
            // Writes X, Y and forces Z = 0.0 (since it's a 2D mesh)
            writeln!(file, "{:.6} {:.6} 0.0", pt.x, pt.y)?;
        }
    }

    println!("Mesh exported to: {}", filename);

    Ok(())
}