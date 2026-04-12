use crate::config;
use ndarray::Array2;
use std::fs::File;
use std::io::{self, Write};


#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

// Initialize mesh directory
pub fn init_mesh_directory(config: &config::Config) {
    let meshname = &config.meshname;
    std::fs::create_dir_all(&format!("job_files/{}/mesh", meshname)).unwrap();
}

// Create grid
pub fn create_grid(longitudinal_points: usize, normal_points: usize) -> Array2<Point> {
    // Inicializa a matriz com tamanho (xi, eta), preenchendo tudo com um Point {x: 0.0, y: 0.0}
    Array2::from_elem((longitudinal_points, normal_points), Point { x: 0.0, y: 0.0 })
}


pub fn save_grid_to_csv(grid: &Array2<Point>, filepath: &str) -> io::Result<()> {
    // Cria (ou sobrescreve) o arquivo no caminho especificado
    let mut file = File::create(filepath)?;

    // Escreve o cabeçalho manualmente
    writeln!(file, "i,j,x,y")?;

    // Itera sobre todos os elementos e seus índices, escrevendo linha por linha
    for ((i, j), point) in grid.indexed_iter() {
        writeln!(file, "{},{},{},{}", i, j, point.x, point.y)?;
    }

    Ok(())
}

/// Thomas algorithm
///
/// Solves a tridiagonal linear system:
///
///     a_j x_{j-1} + b_j x_j + c_j x_{j+1} = d_j
///
/// Used here to solve for all corrections C(i,j) along
/// a vertical line in the Line Gauss-Seidel method.
pub fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {

    let n = d.len();

    let mut c_star = vec![0.0; n];
    let mut d_star = vec![0.0; n];
    let mut x = vec![0.0; n];

    // Forward sweep
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for i in 1..n {
        let denom = b[i] - a[i] * c_star[i - 1];

        if denom.abs() < 1e-14 {
            panic!("Thomas breakdown at i={}", i);
        }

        c_star[i] = if i < n - 1 { c[i] / denom } else { 0.0 };
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / denom;
    }

    // Back substitution
    x[n - 1] = d_star[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    x
}