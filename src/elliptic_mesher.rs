use ndarray::Array2;
use crate::config::Config; // Adjust the import path if necessary
use crate::mesher_utils::Point; // Adjust the import path if necessary
use crate::mesher_utils::thomas_algorithm;
use crate::mesher_utils::periodic_thomas;
use crate::config::EllipticEquation; // Import the new enum for elliptic equation types

/// Entry point for the Elliptic Grid Smoother using an AF1 (Approximate Factorization) scheme.
/// This function solves the Poisson/Laplace equations to produce a smooth, boundary-fitted grid.
/// It uses an Alternating Direction Implicit (ADI) method with multi-frequency acceleration.
pub fn elliptic_smoother(config: &Config, grid: &mut Array2<Point>) {
    println!("Starting elliptic smoothing (ADI)...");

    let ni = config.longitudinal_points - 1; // The index for the overlap (O-mesh branch cut)
    let nj = config.normal_points - 1; // The index for the outer boundary

    // Solver Parameters
    let omega = config.omega; // Over-relaxation factor (typically between 1.0 and 1.8)
    let tolerance = config.conv_criterion;
    let max_iter = config.max_iter;

    // Alpha Sequence (Multi-frequency acceleration)
    // The AF1 scheme uses a sequence of acceleration parameters (alpha) to rapidly 
    // damp out errors across different frequency bandwidths in the grid.
    let alpha_seq = get_alpha_sequence(config, 6); 

    // ==========================================
    // 1. MEMORY PRE-ALLOCATION
    // ==========================================
    // Intermediate Delta (*) arrays for the first step of the ADI scheme
    let mut dx_star = Array2::<f64>::zeros((ni + 1, nj + 1));
    let mut dy_star = Array2::<f64>::zeros((ni + 1, nj + 1));
    
    // Final Delta arrays for the second step of the ADI scheme
    let mut dx_final = Array2::<f64>::zeros((ni + 1, nj + 1));
    let mut dy_final = Array2::<f64>::zeros((ni + 1, nj + 1));

    // Vectors for Periodic Thomas Algorithm (Xi direction / wrap-around) (Size ni)
    let mut a_xi = vec![0.0; ni];
    let mut b_xi = vec![0.0; ni];
    let mut c_xi = vec![0.0; ni];
    let mut d_x_xi = vec![0.0; ni];
    let mut d_y_xi = vec![0.0; ni];

    // Vectors for Standard Thomas Algorithm (Eta direction / wall-to-farfield) 
    // Size is nj - 1, since j=0 (wall) and j=nj (farfield) are fixed Dirichlet boundaries
    let n_j_internal = nj - 1;
    let mut a_eta = vec![0.0; n_j_internal];
    let mut b_eta = vec![0.0; n_j_internal];
    let mut c_eta = vec![0.0; n_j_internal];
    let mut d_x_eta = vec![0.0; n_j_internal];
    let mut d_y_eta = vec![0.0; n_j_internal];

    // ==========================================
    // 2. GLOBAL CONVERGENCE LOOP
    // ==========================================
    let mut iter = 0;
    let mut max_error = 1.0;

    while max_error > tolerance && iter < max_iter {
        max_error = 0.0;

        // Step 1: Handle Boundary Control Terms (Sorenson & Steger technique)
        // If Poisson is selected, compute the P0 and Q0 source terms at the wall to 
        // strictly enforce grid orthogonality and initial normal spacing (ds_wall).
        let (p0, q0) = match config.elliptic_equation {
            EllipticEquation::Poisson { ds_wall, .. } => compute_boundary_controls(grid, ni, ds_wall, iter),
            EllipticEquation::Laplace => (vec![0.0; ni], vec![0.0; ni]),
        };


        for &alpha in &alpha_seq {
            
            // --- STEP 1: XI SWEEP (Longitudinal / Periodic) ---
            // Solves the first part of the factored operator: (alpha - A_xi) * Delta* = alpha * omega * Residual
            for j in 1..nj {
                for i in 0..ni {
                    // Calculate metric tensor components and PDE residuals using central finite differences
                    let (a_ij, _, _, res_x, res_y) = calculate_metrics_and_residual(config, grid, i, j, ni, &p0, &q0);

                    // Tridiagonal matrix coefficients for the Xi Sweep
                    a_xi[i] = -a_ij;
                    b_xi[i] = 2.0 * a_ij + alpha; // Alpha acts as a diagonal dominant stabilizer
                    c_xi[i] = -a_ij;
                    
                    // RHS vectors
                    d_x_xi[i] = alpha * omega * res_x;
                    d_y_xi[i] = alpha * omega * res_y;
                }

                // Solve the periodic tridiagonal system for the current radial ring 'j'
                let out_x_xi = periodic_thomas(&a_xi, &b_xi, &c_xi, &d_x_xi);
                let out_y_xi = periodic_thomas(&a_xi, &b_xi, &c_xi, &d_y_xi);

                // Save the intermediate correction result (\Delta^*)
                for i in 0..ni {
                    dx_star[[i, j]] = out_x_xi[i];
                    dy_star[[i, j]] = out_y_xi[i];
                }
            }

            // --- STEP 2: ETA SWEEP (Normal / Dirichlet Fixed Boundaries) ---
            // Solves the second part of the factored operator: (alpha - A_eta) * Delta = alpha * Delta*
            for i in 0..ni {
                for j_int in 0..n_j_internal {
                    let j = j_int + 1; // Maps the vector index (0) to the matrix index (1)
                    
                    let (_, _, c_ij, _, _) = calculate_metrics_and_residual(config, grid, i, j, ni, &p0, &q0 );

                    // Tridiagonal matrix coefficients for the Eta Sweep
                    a_eta[j_int] = -c_ij;
                    b_eta[j_int] = 2.0 * c_ij + alpha;
                    c_eta[j_int] = -c_ij;
                    
                    // The right-hand side is the intermediate matrix scaled by alpha
                    d_x_eta[j_int] = alpha * dx_star[[i, j]];
                    d_y_eta[j_int] = alpha * dy_star[[i, j]];
                }

                // Solve the standard tridiagonal system along the normal line 'i' (airfoil to far-field)
                let out_x_eta = thomas_algorithm(&a_eta, &b_eta, &c_eta, &d_x_eta);
                let out_y_eta = thomas_algorithm(&a_eta, &b_eta, &c_eta, &d_y_eta);

                // Save the final correction result and calculate the maximum residual
                for j_int in 0..n_j_internal {
                    let j = j_int + 1;
                    dx_final[[i, j]] = out_x_eta[j_int];
                    dy_final[[i, j]] = out_y_eta[j_int];

                    // Track maximum displacement for convergence checking
                    let current_error = dx_final[[i, j]].abs().max(dy_final[[i, j]].abs());
                    if current_error > max_error {
                        max_error = current_error;
                    }
                }
            }

            // --- 3. MESH UPDATE ---
            // Apply the final delta corrections to the physical grid coordinates
            for j in 1..nj {
                for i in 0..ni {
                    grid[[i, j]].x += dx_final[[i, j]];
                    grid[[i, j]].y += dy_final[[i, j]];
                }
                // Ensure perfect periodic closure of the O-mesh at the branch cut (wake line)
                grid[[ni, j]].x = grid[[0, j]].x;
                grid[[ni, j]].y = grid[[0, j]].y;
            }
        } // End of the Alpha Sequence loop

        iter += 1;
        if iter % 100 == 0 || iter == 1 {
            println!("Iteration: {}/{}, Max Residual: {:.2e}", iter, max_iter, max_error);
        }
        
        // Divergence check (if the solution explodes, typically due to aggressive source terms or stretching)
        if max_error > 1e5 {
            panic!("Solution diverged at iteration: {}, Max Residual: {:.2e} \n", iter, max_error);
        }

    }

    println!("Smoothing completed in {} global iterations.", iter);
}

/// Generates a geometric sequence of acceleration parameters (alpha).
/// This multi-grid-like approach cycles through small and large alphas to damp 
/// high-frequency (local) and low-frequency (global) numerical errors efficiently.
fn get_alpha_sequence(config: &Config, m_total: usize) -> Vec<f64> {
    let alpha_l = config.alpha_l; // Min (stability / global error damping)
    let alpha_h = config.alpha_h; // Max (acceleration / local error damping)

    if m_total <= 1 {
        return vec![alpha_l];
    }

    let mut alpha_seq = vec![0.0; m_total];
    let denominator = m_total as f64 - 1.0;

    for m in 0..m_total {
        let fraction = (m as f64) / denominator;
        
        let alpha_m = alpha_h * (alpha_l / alpha_h).powf(fraction);
        alpha_seq[m] = alpha_m;
    }

    // Applying from highest to lowest is usually standard practice for AF schemes
    alpha_seq.reverse();
    alpha_seq
        
}

/// Helper function to compute the covariant metric tensor elements and the PDE residuals.
/// Returns a tuple: (a_ij, b_ij, c_ij, res_x, res_y)
fn calculate_metrics_and_residual(
                                  config: &Config, grid: &ndarray::Array2<Point>, i: usize, j: usize, 
                                  ni: usize, p0: &[f64], q0: &[f64]
                                  ) -> (f64, f64, f64, f64, f64) {    
    
    // Periodic treatment for neighbors in 'i'
    let i_m1 = if i == 0 { ni - 1 } else { i - 1 };
    let i_p1 = if i == ni - 1 { 0 } else { i + 1 };

    // Allocate variables for Poisson source terms
    let mut p_ij = 0.0;
    let mut q_ij = 0.0;

    // Neighbors in 'j' (The caller loop boundaries ensure we never access j < 0 or j > nj)
    let j_m1 = j - 1;
    let j_p1 = j + 1;

    // Coordinates of central node and cross-neighbors
    let x_ij = grid[[i, j]].x;     let y_ij = grid[[i, j]].y;
    let x_ip1 = grid[[i_p1, j]].x; let y_ip1 = grid[[i_p1, j]].y;
    let x_im1 = grid[[i_m1, j]].x; let y_im1 = grid[[i_m1, j]].y;
    let x_jp1 = grid[[i, j_p1]].x; let y_jp1 = grid[[i, j_p1]].y;
    let x_jm1 = grid[[i, j_m1]].x; let y_jm1 = grid[[i, j_m1]].y;

    // First Derivatives (Central Differences in the computational domain)
    let x_xi  = 0.5 * (x_ip1 - x_im1);
    let y_xi  = 0.5 * (y_ip1 - y_im1);
    let x_eta = 0.5 * (x_jp1 - x_jm1);
    let y_eta = 0.5 * (y_jp1 - y_jm1);

    // Metric Tensor Coefficients (often denoted g22, g12, g11 or Alpha, Beta, Gamma in literature)
    // a_ij controls diffusion along xi, c_ij controls diffusion along eta
    let a_ij = x_eta * x_eta + y_eta * y_eta;
    let b_ij = x_xi * x_eta + y_xi * y_eta; 
    let c_ij = x_xi * x_xi + y_xi * y_xi;
    
    // Second Derivatives
    let x_xixi = x_ip1 - 2.0 * x_ij + x_im1;
    let y_xixi = y_ip1 - 2.0 * y_ij + y_im1;
    let x_etaeta = x_jp1 - 2.0 * x_ij + x_jm1;
    let y_etaeta = y_jp1 - 2.0 * y_ij + y_jm1;
    
    // Cross Derivative (x_xieta) - Grabs the 4 diagonal corners of the stencil
    let x_ip1_jp1 = grid[[i_p1, j_p1]].x; let y_ip1_jp1 = grid[[i_p1, j_p1]].y;
    let x_im1_jm1 = grid[[i_m1, j_m1]].x; let y_im1_jm1 = grid[[i_m1, j_m1]].y;
    let x_im1_jp1 = grid[[i_m1, j_p1]].x; let y_im1_jp1 = grid[[i_m1, j_p1]].y;
    let x_ip1_jm1 = grid[[i_p1, j_m1]].x; let y_ip1_jm1 = grid[[i_p1, j_m1]].y;
    
    let x_xieta = 0.25 * (x_ip1_jp1 - x_ip1_jm1 - x_im1_jp1 + x_im1_jm1);
    let y_xieta = 0.25 * (y_ip1_jp1 - y_ip1_jm1 - y_im1_jp1 + y_im1_jm1);
    
    // Control functions P and Q terms
    let d_ij = (x_xi * y_eta - x_eta * y_xi).powi(2); // Squared Jacobian determinant

    // Calculate decaying exponential source terms based on wall boundary values
    // This allows the grid control (orthogonality/spacing) to be extremely strong at the wall,
    // but blend smoothly back to Laplace (uncontrolled smoothing) near the far field.
    if let EllipticEquation::Poisson { a_decay, b_decay, .. } = config.elliptic_equation {
        let dist = j as f64; // Distance in computational space
        p_ij = p0[i] * (-a_decay * dist).exp();
        q_ij = q0[i] * (-b_decay * dist).exp();
    }

    // Elliptic Residuals 
    // This is the discretized evaluation of: a*r_xixi - 2b*r_xieta + c*r_etaeta + J^2*(P*r_xi + Q*r_eta)
    let res_x = a_ij * x_xixi - 2.0 * b_ij * x_xieta + c_ij * x_etaeta + d_ij * (p_ij * x_xi + q_ij * x_eta);
    let res_y = a_ij * y_xixi - 2.0 * b_ij * y_xieta + c_ij * y_etaeta + d_ij * (p_ij * y_xi + q_ij * y_eta);

    (a_ij, b_ij, c_ij, res_x, res_y)
}

/// Computes the Steger-Sorenson source terms (P0 and Q0) dynamically at the inner boundary (j=0).
/// These terms force the generating Poisson equation to push lines into an orthogonal state 
/// relative to the wall, while matching a strictly defined normal spacing distance (ds_wall).
fn compute_boundary_controls(
    grid: &ndarray::Array2<Point>, 
    ni: usize, 
    ds_wall: f64,
    iter: usize
) -> (Vec<f64>, Vec<f64>) {
    let mut p0 = vec![0.0; ni];
    let mut q0 = vec![0.0; ni];

    for i in 0..ni {
        let i_m1 = if i == 0 { ni - 1 } else { i - 1 };
        let i_p1 = if i == ni - 1 { 0 } else { i + 1 };

        // 1. Surface tangent derivatives (j=0)
        let x_xi = 0.5 * (grid[[i_p1, 0]].x - grid[[i_m1, 0]].x);
        let y_xi = 0.5 * (grid[[i_p1, 0]].y - grid[[i_m1, 0]].y);
        let s1 = (x_xi.powi(2) + y_xi.powi(2)).sqrt();

        // Prevent division by zero if nodes are collapsed
        if s1 < 1e-12 { continue; }

        // 2. Desired normal derivatives (Strict orthogonality condition + ds_wall spacing)
        let x_eta_des = -ds_wall * (y_xi / s1);
        let y_eta_des =  ds_wall * (x_xi / s1);

        // 3. Metric coefficients at the wall forced by the desired target derivatives
        let a_f = x_eta_des.powi(2) + y_eta_des.powi(2);
        let b_f = 0.0; // Enforced perfect orthogonality (dot product of xi and eta vectors = 0)
        let c_f = s1.powi(2);
        
        let j_inv = x_xi * y_eta_des - x_eta_des * y_xi;
        let d_f = j_inv.powi(2);

        // 4. Second derivatives (Stencil isolated at the boundary)
        let x_xi_xi = grid[[i_p1, 0]].x - 2.0 * grid[[i, 0]].x + grid[[i_m1, 0]].x;
        let y_xi_xi = grid[[i_p1, 0]].y - 2.0 * grid[[i, 0]].y + grid[[i_m1, 0]].y;

        // One-sided second derivatives for eta at j=0 (using interior points j=0, 1, 2)
        let x_eta_eta = grid[[i, 2]].x - 2.0 * grid[[i, 1]].x + grid[[i, 0]].x;
        let y_eta_eta = grid[[i, 2]].y - 2.0 * grid[[i, 1]].y + grid[[i, 0]].y;

        // Cross derivatives: central difference of xi-derivatives between j=1 and j=0
        let x_xi_j1 = 0.5 * (grid[[i_p1, 1]].x - grid[[i_m1, 1]].x);
        let y_xi_j1 = 0.5 * (grid[[i_p1, 1]].y - grid[[i_m1, 1]].y);
        let x_xi_eta = x_xi_j1 - x_xi;
        let y_xi_eta = y_xi_j1 - y_xi;

        // 5. Solve the 2x2 system via Cramer's Rule to find the required P0 and Q0 
        // that force the residual of the PDE to zero under these geometric conditions.
        let rhs_x = -(a_f * x_xi_xi - 2.0 * b_f * x_xi_eta + c_f * x_eta_eta) / d_f;
        let rhs_y = -(a_f * y_xi_xi - 2.0 * b_f * y_xi_eta + c_f * y_eta_eta) / d_f;

        // Relaxation factor: Because Steger-Sorenson source terms are highly non-linear, 
        // applying them fully at iter 1 will cause the mesh to explode. We ramp them up 
        // smoothly over the first 1000 iterations.
        let p_relax = 0.005 * (iter as f64 / 1000.0); 
        let q_relax = 0.001 * (iter as f64 / 1000.0); 

        p0[i] = (y_eta_des * rhs_x - x_eta_des * rhs_y) / j_inv * p_relax;
        q0[i] = (-y_xi * rhs_x + x_xi * rhs_y) / j_inv * q_relax;
    }
    (p0, q0)
}