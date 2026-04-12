use crate::config::Config;
use ndarray::Array2;
use crate::mesher_utils::Point;
use crate::mesher_utils::thomas;


pub fn generate_parabolic_mesh(config: &Config, grid: &mut Array2<Point>) {
    println!("Generating parabolic mesh...");
    
    space_march(config, grid);

}

fn space_march(config: &Config, grid: &mut Array2<Point>) {
    let longitudinal_points = config.longitudinal_points;
    let normal_points = config.normal_points;

    for j in 1..normal_points-1 {

        // Alocate vectors 
        let mut ref_grid_jp1 = vec![Point { x: 0.0, y: 0.0 }; longitudinal_points];
        let mut a = vec![0.0; longitudinal_points];
        let mut b = vec![0.0; longitudinal_points];
        let mut c = vec![0.0; longitudinal_points];
        let mut d_x = vec![0.0; longitudinal_points];
        let mut d_y = vec![0.0; longitudinal_points];

        // Calculate the local reference grid point for each i at j+1
        let s_j = streching_factor(config, j);
        let s_jm1 = streching_factor(config, j-1);
        let s_jnj = streching_factor(config, config.normal_points - 1);
        let ds = (s_j - s_jm1) / (s_jnj - s_jm1);

        // Calculate the local reference grid points for all i at j+1
        for i in 1..longitudinal_points-1 {
            ref_grid_jp1[i] = local_reference_grid(config, grid, i, j+1, ds);
        }

        for i in 1..longitudinal_points-1 {
            // xi derivatives
            let x_xi = 1.0 / 2.0 * (grid[[i+1, j-1]].x - grid[[i-1, j-1]].x);
            let y_xi = 1.0 / 2.0 * (grid[[i+1, j-1]].y - grid[[i-1, j-1]].y);

            // eta derivatives
            let x_eta = ref_grid_jp1[i].x - grid[[i, j]].x;
            let y_eta = ref_grid_jp1[i].y - grid[[i, j]].y;

            let A = x_eta.powi(2) + y_eta.powi(2);
            let B = x_xi * x_eta + y_xi * y_eta;
            let C = x_xi.powi(2) + y_xi.powi(2);

            // Coefficients for the tridiagonal system
            a[i] = 2.0 * A;
            b[i] = -4.0 * (A + C);
            c[i] = 2.0 * A;
            
            // LHS of the tridiagonal system
            d_x[i] = B * (ref_grid_jp1[i+1].x - ref_grid_jp1[i-1].x - grid[[i+1, j-1]].x + grid[[i-1, j-1]].x) -
                        2.0 * C * (ref_grid_jp1[i].x - grid[[i, j-1]].x); 
            d_y[i] = B * (ref_grid_jp1[i+1].y - ref_grid_jp1[i-1].y - grid[[i+1, j-1]].y + grid[[i-1, j-1]].y) -
                        2.0 * C * (ref_grid_jp1[i].y - grid[[i, j-1]].y); 
        }

        // Solve tridiagonal system to obtain corrections C(i,j)
        let x = thomas(&a, &b, &c, &d_x);
        let y = thomas(&a, &b, &c, &d_y);

        for i in 1..longitudinal_points-1 {
            grid[[i, j]].x = x[i];
            grid[[i, j]].y = y[i];
        }
    }

}

fn local_reference_grid(config: &Config, grid: &mut Array2<Point>, i: usize, j: usize, ds: f64) -> Point{
    // outter boundary index
    let nj = config.normal_points - 1;

    // xi derivatives
    let x_xi = 1.0 / 2.0 * (grid[[i+1, j-1]].x - grid[[i-1, j-1]].x);
    let y_xi = 1.0 / 2.0 * (grid[[i+1, j-1]].y - grid[[i-1, j-1]].y);
    
    let S_eta = calc_S_eta(config, grid, i, j, ds);

    // eta derivatives
    let x_eta = - (y_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();
    let y_eta = (x_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();

    // grid point coord
    let x_o = grid[[i, j-1]].x + x_eta;
    let y_o = grid[[i, j-1]].y + y_eta;

    // interpolated points
    let x_int = grid[[i, j-1]].x + ds * (grid[[i, nj]].x - grid[[i, j-1]].x);
    let y_int = grid[[i, j-1]].y + ds * (grid[[i, nj]].y - grid[[i, j-1]].y);

    // Normalized distance from the inner boundary (j=0) to the outer boundary (j=nj)
    let eps = (j as f64 - 1.0) / (nj as f64 - 1.0);

    // local reference grid point
    let x_ref = eps * x_int + (1.0 - eps) * x_o;
    let y_ref = eps * y_int + (1.0 - eps) * y_o;

    Point { x: x_ref, y: y_ref }
}

fn streching_factor(config: &Config, j: usize) -> f64 {
    // outter boundary index
    let nj = config.normal_points - 1;
    // Mesh parameters
    let P = config.P;
    let Q = config.Q;

    // Converte j e nj para f64 ANTES da divisão para garantir a precisão decimal.
    // Assumindo que j varia de 0 até (nj - 1).
    let eta_norm = (j as f64) / (nj as f64);

    // Aplicação da fórmula usando o método .tanh()
    let s = P * eta_norm + (1.0 - P) * (1.0 - ((Q * (1.0 - eta_norm)).tanh() / Q.tanh()));

    s
}

fn calc_S_eta(config: &Config, grid: &mut Array2<Point>,i: usize, j: usize, ds: f64) -> f64 {
    // outter boundary index
    let nj = config.normal_points - 1;

    let dx = grid[[i, j-1]].x - grid[[i, nj]].x;
    let dy = grid[[i, j-1]].y - grid[[i, nj]].y;
    let R = (dx.powi(2) + dy.powi(2)).sqrt();

    // Aplicação da fórmula usando o método .tanh()
    let S_eta = ds * R;

    S_eta
}