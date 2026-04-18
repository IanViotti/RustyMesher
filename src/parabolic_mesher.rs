use std::vec;

use crate::config::Config;
use ndarray::Array2;
use crate::mesher_utils::Point;
use crate::mesher_utils::periodic_thomas;


pub fn generate_parabolic_mesh(config: &Config, grid: &mut Array2<Point>) {
    println!("Generating parabolic mesh...");
    
    space_march(config, grid);

}

fn space_march(config: &Config, grid: &mut Array2<Point>) {
    let ni = config.longitudinal_points - 1;
    let nj = config.normal_points - 1;

    let s_bar = streching_factor(config);

    for j in 1..nj{
        println!("Processing layer j: {}", j);

        // Alocate vectors 
        //let mut ref_grid_jp1: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];
        //let mut ref_grid_j: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];
        let mut a = vec![0.0; ni];
        let mut b = vec![0.0; ni];
        let mut c = vec![0.0; ni];
        let mut d_x = vec![0.0; ni];
        let mut d_y = vec![0.0; ni];

        // Grid spacing control in eta direction
        let ds = (s_bar[j] - s_bar[j-1]) / (s_bar[nj-1] - s_bar[j-1]);

        // Calculate the local reference grid points for all i at j+1
        let (ref_grid_j, ref_grid_jp1) = local_reference_grid(config, grid, j, ds);

        println!("j: local ref grid {} -{:?}", j, ref_grid_j);

        for i in 0..ni{

            // i indexes for periodic boundary conditions
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

            // xi derivatives
            let x_xi = 1.0 / 2.0 * (ref_grid_jp1[ip1].x - ref_grid_jp1[im1].x);
            let y_xi = 1.0 / 2.0 * (ref_grid_jp1[ip1].y - ref_grid_jp1[im1].y);

            // eta derivatives
            let x_eta = ref_grid_jp1[i].x - ref_grid_j[i].x;
            let y_eta = ref_grid_jp1[i].y - ref_grid_j[i].y;

            let A = x_eta.powi(2) + y_eta.powi(2);
            let B = x_xi * x_eta + y_xi * y_eta;
            let C = x_xi.powi(2) + y_xi.powi(2);

            // Coefficients for the tridiagonal system
            a[i] = 2.0 * A;
            b[i] = -4.0 * (A + C);
            c[i] = 2.0 * A;
            
            // LHS of the tridiagonal system
            d_x[i] = B * (ref_grid_jp1[ip1].x - ref_grid_jp1[im1].x - grid[[ip1, j-1]].x + grid[[im1, j-1]].x) -
                        2.0 * C * (ref_grid_jp1[i].x - grid[[i, j-1]].x); 
            d_y[i] = B * (ref_grid_jp1[ip1].y - ref_grid_jp1[im1].y - grid[[ip1, j-1]].y + grid[[im1, j-1]].y) -
                        2.0 * C * (ref_grid_jp1[i].y - grid[[i, j-1]].y); 
        }

        // Solve tridiagonal system to obtain corrections C(i,j)
        let x = periodic_thomas(&a, &b, &c, &d_x);
        let y = periodic_thomas(&a, &b, &c, &d_y);

        for i in 0..ni{
            grid[[i, j]].x = x[i]; 
            grid[[i, j]].y = y[i];

            //grid[[i, j]] = ref_grid_jp1[i]; 

            println!("grid @ i:{}, j:{} = {:?}", i, j, grid[[i, j]])
        }
        // pontos sobrepostos para garantir periodicidade
        grid[[ni, j]] = grid[[0, j]];

    }
}


// Calculates the local reference grid @ j and j+1
// Returns a tuple with vectors for j and j+1
fn local_reference_grid(config: &Config, grid: &mut Array2<Point>,j: usize, ds: f64) -> (Vec<Point>, Vec<Point>) {

    // outter boundary index
    let nj = config.normal_points - 1;
    let ni = config.longitudinal_points - 1;

    // Alocate vector 
    let mut ref_grid_j: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];
    let mut ref_grid_jp1: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }; ni];

    for i in 0..ni-1{
    
        // i indexes for periodic boundary conditions
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

        // xi derivatives
        let x_xi = 1.0 / 2.0 * (grid[[ip1, j-1]].x - grid[[im1, j-1]].x);
        let y_xi = 1.0 / 2.0 * (grid[[ip1, j-1]].y - grid[[im1, j-1]].y);

        let dx = grid[[i, nj]].x - grid[[i, j-1]].x;
        let dy = grid[[i, nj]].y - grid[[i, j-1]].y;
        let S_eta = calc_S_eta(dx, dy, ds);

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
        ref_grid_j[i].x = eps * x_int + (1.0 - eps) * x_o;
        ref_grid_j[i].y = eps * y_int + (1.0 - eps) * y_o;

    }


    // === Local Reference Grid @ j + 1 ===

    for i in 0..ni-1{

        // i indexes for periodic boundary conditions
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

        // xi derivatives
        let x_xi = 1.0 / 2.0 * (ref_grid_j[ip1].x - ref_grid_j[im1].x);
        let y_xi = 1.0 / 2.0 * (ref_grid_j[ip1].y - ref_grid_j[im1].y);

        let dx = grid[[i, nj]].x - ref_grid_j[i].x;
        let dy = grid[[i, nj]].y - ref_grid_j[i].y;
        let S_eta = calc_S_eta(dx, dy, ds);

        // eta derivatives
        let x_eta = - (y_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();
        let y_eta = (x_xi * S_eta) / (x_xi.powi(2) + y_xi.powi(2)).sqrt();

        // grid point coord
        let x_o = ref_grid_j[i].x + x_eta;
        let y_o = ref_grid_j[i].y + y_eta;

        // interpolated points
        let x_int = ref_grid_j[i].x + ds * (grid[[i, nj]].x - ref_grid_j[i].x);
        let y_int = ref_grid_j[i].y + ds * (grid[[i, nj]].y - ref_grid_j[i].y);

        // Normalized distance from the inner boundary (j=0) to the outer boundary (j=nj)
        let eps = (j as f64 - 1.0) / (nj as f64 - 1.0);

        // local reference grid point
        ref_grid_jp1[i].x = eps * x_int + (1.0 - eps) * x_o;
        ref_grid_jp1[i].y = eps * y_int + (1.0 - eps) * y_o;
    }

    (ref_grid_j, ref_grid_jp1)

}

fn streching_factor(config: &Config) -> Vec<f64> {
    // outter boundary index
    let nj = config.normal_points - 1;
    // Mesh parameters
    let q = config.q;

    let mut s_bar = vec![0.0; nj];

    let mut s_bar_jm1: f64;
    let mut s_bar_j: f64;

    for j in 0..nj-1{
        if j == 0 {
            s_bar_jm1 = 0.0;
            s_bar_j = q;
        }
        else {
            s_bar_jm1 = s_bar[j-1];
            s_bar_j = s_bar[j];
        }
        s_bar[j+1] = s_bar_j + q * (s_bar_j - s_bar_jm1)
    }

    s_bar
}

fn calc_S_eta(dx: f64, dy: f64, ds: f64) -> f64 {

    let R = (dx.powi(2) + dy.powi(2)).sqrt();

    // Aplicação da fórmula usando o método .tanh()
    let S_eta = ds * R;

    S_eta
}