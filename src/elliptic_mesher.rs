use ndarray::Array2;
use crate::config::Config; // Adjust the import path if necessary
use crate::mesher_utils::Point; // Adjust the import path if necessary
use crate::mesher_utils::thomas_algorithm;
use crate::mesher_utils::periodic_thomas;

pub fn elliptic_smoother(config: &Config, grid: &mut Array2<Point>) {
    println!("Starting elliptic smoothing (ADI)...");

    let ni = config.longitudinal_points - 1; // The index for the overlap (O-mesh branch cut)
    let nj = config.normal_points - 1; // The index for the outer boundary

    // Solver Parameters
    let omega = config.omega; // Relaxation factor (typically between 1.0 and 1.8)
    let tolerance = config.conv_criterion;
    let max_iter = config.max_iter;

    // Alpha Sequence (Multi-frequency acceleration)
    // Typical values based on literature for O-meshes:
    let alpha_seq = get_alpha_sequence(config, 6); 

    // ==========================================
    // 1. MEMORY PRE-ALLOCATION
    // ==========================================
    let mut dx_star = Array2::<f64>::zeros((ni + 1, nj + 1));
    let mut dy_star = Array2::<f64>::zeros((ni + 1, nj + 1));
    let mut dx_final = Array2::<f64>::zeros((ni + 1, nj + 1));
    let mut dy_final = Array2::<f64>::zeros((ni + 1, nj + 1));

    // Vectors for Periodic Thomas (Size ni)
    let mut a_xi = vec![0.0; ni];
    let mut b_xi = vec![0.0; ni];
    let mut c_xi = vec![0.0; ni];
    let mut d_x_xi = vec![0.0; ni];
    let mut d_y_xi = vec![0.0; ni];

    // Vectors for Standard Thomas (Size nj - 1, since j=0 and j=nj are fixed boundaries)
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

        for &alpha in &alpha_seq {
            
            // --- STEP 1: XI SWEEP (Longitudinal / Periodic) ---
            for j in 1..nj {
                for i in 0..ni {
                    // Calculate metrics using central finite differences on the current grid
                    let (a_ij, _, _, res_x, res_y) = calculate_metrics_and_residual(grid, i, j, ni);

                    // SymPy coefficients for the Xi Sweep
                    a_xi[i] = -a_ij;
                    b_xi[i] = 2.0 * a_ij + alpha;
                    c_xi[i] = -a_ij;
                    d_x_xi[i] = alpha * omega * res_x;
                    d_y_xi[i] = alpha * omega * res_y;
                }

                // Solve the periodic tridiagonal system for the current ring 'j'
                let out_x_xi = periodic_thomas(&a_xi, &b_xi, &c_xi, &d_x_xi);
                let out_y_xi = periodic_thomas(&a_xi, &b_xi, &c_xi, &d_y_xi);

                // Save the intermediate result
                for i in 0..ni {
                    dx_star[[i, j]] = out_x_xi[i];
                    dy_star[[i, j]] = out_y_xi[i];
                }
            }

            // --- STEP 2: ETA SWEEP (Normal / Dirichlet Fixed Boundaries) ---
            for i in 0..ni {
                for j_int in 0..n_j_internal {
                    let j = j_int + 1; // Maps the vector index (0) to the matrix index (1)
                    
                    let (_, _, c_ij, _, _) = calculate_metrics_and_residual(grid, i, j, ni);

                    // SymPy coefficients for the Eta Sweep
                    a_eta[j_int] = -c_ij;
                    b_eta[j_int] = 2.0 * c_ij + alpha;
                    c_eta[j_int] = -c_ij;
                    
                    // The right-hand side is the intermediate matrix multiplied by alpha!
                    d_x_eta[j_int] = alpha * dx_star[[i, j]];
                    d_y_eta[j_int] = alpha * dy_star[[i, j]];
                }

                // Solve the standard tridiagonal system along the radial line 'i'
                let out_x_eta = thomas_algorithm(&a_eta, &b_eta, &c_eta, &d_x_eta);
                let out_y_eta = thomas_algorithm(&a_eta, &b_eta, &c_eta, &d_y_eta);

                // Save the final result and calculate the maximum error
                for j_int in 0..n_j_internal {
                    let j = j_int + 1;
                    dx_final[[i, j]] = out_x_eta[j_int];
                    dy_final[[i, j]] = out_y_eta[j_int];

                    let current_error = dx_final[[i, j]].abs().max(dy_final[[i, j]].abs());
                    if current_error > max_error {
                        max_error = current_error;
                    }
                }
            }

            // --- 3. MESH UPDATE ---
            for j in 1..nj {
                for i in 0..ni {
                    grid[[i, j]].x += dx_final[[i, j]];
                    grid[[i, j]].y += dy_final[[i, j]];
                }
                // Ensure perfect periodic closure of the O-mesh (branch cut)
                grid[[ni, j]].x = grid[[0, j]].x;
                grid[[ni, j]].y = grid[[0, j]].y;
            }
        } // End of the Alpha Sequence loop

        iter += 1;
        if iter % 100 == 0 {
            println!("Iteration: {}/{}, Max Residual: {:.2e}", iter, max_iter, max_error);
        }

        // Divergence check if the solution overflows
        if max_error > 1e5 {
            panic!("Solution diverged at iteration: {}, Max Residual: {:.2e} \n", iter, max_error);
        }

    }

    println!("Smoothing completed in {} global iterations.", iter);
}

fn get_alpha_sequence(config: &Config, m_total: usize) -> Vec<f64> {
    // Geometric progression for alpha values
    let alpha_l = config.alpha_l;  // Min (stability)
    let alpha_h = config.alpha_h; // max (acceleration)

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

    alpha_seq.reverse();
    alpha_seq
        
}

/// Helper function to compute the transformation metrics and the PDE residual
fn calculate_metrics_and_residual(grid: &Array2<Point>, i: usize, j: usize, ni: usize) -> (f64, f64, f64, f64, f64) {
    // Periodic treatment for neighbors in 'i'
    let i_m1 = if i == 0 { ni - 1 } else { i - 1 };
    let i_p1 = if i == ni - 1 { 0 } else { i + 1 };

    // Neighbors in 'j' (The loop boundaries ensure we never access out-of-bounds indices)
    let j_m1 = j - 1;
    let j_p1 = j + 1;

    // Coordinates
    let x_ij = grid[[i, j]].x;     let y_ij = grid[[i, j]].y;
    let x_ip1 = grid[[i_p1, j]].x; let y_ip1 = grid[[i_p1, j]].y;
    let x_im1 = grid[[i_m1, j]].x; let y_im1 = grid[[i_m1, j]].y;
    let x_jp1 = grid[[i, j_p1]].x; let y_jp1 = grid[[i, j_p1]].y;
    let x_jm1 = grid[[i, j_m1]].x; let y_jm1 = grid[[i, j_m1]].y;

    // First Derivatives (Central Differences)
    let x_xi  = 0.5 * (x_ip1 - x_im1);
    let y_xi  = 0.5 * (y_ip1 - y_im1);
    let x_eta = 0.5 * (x_jp1 - x_jm1);
    let y_eta = 0.5 * (y_jp1 - y_jm1);

    // Metric Coefficients (Alpha, Beta, Gamma from the literature)
    let a_ij = x_eta * x_eta + y_eta * y_eta;
    let b_ij = x_xi * x_eta + y_xi * y_eta; 
    let c_ij = x_xi * x_xi + y_xi * y_xi;

    // Second Derivatives
    let x_xixi = x_ip1 - 2.0 * x_ij + x_im1;
    let y_xixi = y_ip1 - 2.0 * y_ij + y_im1;
    let x_etaeta = x_jp1 - 2.0 * x_ij + x_jm1;
    let y_etaeta = y_jp1 - 2.0 * y_ij + y_jm1;

    // Cross Derivative (x_xieta) - Grabs the 4 diagonal corners
    let x_ip1_jp1 = grid[[i_p1, j_p1]].x; let y_ip1_jp1 = grid[[i_p1, j_p1]].y;
    let x_im1_jm1 = grid[[i_m1, j_m1]].x; let y_im1_jm1 = grid[[i_m1, j_m1]].y;
    let x_im1_jp1 = grid[[i_m1, j_p1]].x; let y_im1_jp1 = grid[[i_m1, j_p1]].y;
    let x_ip1_jm1 = grid[[i_p1, j_m1]].x; let y_ip1_jm1 = grid[[i_p1, j_m1]].y;

    let x_xieta = 0.25 * (x_ip1_jp1 - x_ip1_jm1 - x_im1_jp1 + x_im1_jm1);
    let y_xieta = 0.25 * (y_ip1_jp1 - y_ip1_jm1 - y_im1_jp1 + y_im1_jm1);

    // Residuals (Notice the subtraction of the B term since it is defined as positive)
    let res_x = a_ij * x_xixi - 2.0 * b_ij * x_xieta + c_ij * x_etaeta;
    let res_y = a_ij * y_xixi - 2.0 * b_ij * y_xieta + c_ij * y_etaeta;

    (a_ij, b_ij, c_ij, res_x, res_y)
}