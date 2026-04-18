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
pub fn thomas_algorithm(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let mut c_star = vec![0.0; n];
    let mut d_star = vec![0.0; n];
    let mut x = vec![0.0; n];

    // Passo 1: Eliminação (Forward sweep) - Construindo L e U implicitamente
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for i in 1..n {
        let m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        if i < n - 1 {
            c_star[i] = c[i] * m;
        }
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    // Passo 2: Substituição regressiva (Backward substitution)
    x[n - 1] = d_star[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }

    x
}

/// Resolve um sistema tridiagonal cíclico/periódico via Sherman-Morrison
pub fn periodic_thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    
    // Define o gama (usualmente -b[0] para garantir estabilidade numérica)
    let gamma = -b[0]; 
    
    // Cria a diagonal principal modificada (Matriz T)
    let mut bb = b.to_vec();
    bb[0] = b[0] - gamma;
    bb[n - 1] = b[n - 1] - a[0] * c[n - 1] / gamma;

    // Monta o vetor u (apenas as pontas possuem valores)
    let mut u = vec![0.0; n];
    u[0] = gamma;
    u[n - 1] = a[0]; 

    // O vetor v é modelado implicitamente para economizar alocação:
    // v = [1.0, 0.0, ..., 0.0, c[n-1]/gamma]

    // Resolução Dupla:
    // 1. Resolve Ty = d usando a diagonal principal modificada
    let y: Vec<f64> = thomas_algorithm(a, &bb, c, d);
    
    // 2. Resolve Tq = u
    let q = thomas_algorithm(a, &bb, c, &u);

    // 3. Reconstrução final da solução combinando os resultados
    let v_n_minus_1 = c[n - 1] / gamma;
    
    // Produtos escalares (v dot y) e (v dot q)
    let v_dot_y = y[0] + v_n_minus_1 * y[n - 1];
    let v_dot_q = q[0] + v_n_minus_1 * q[n - 1];

    let factor = v_dot_y / (1.0 + v_dot_q);

    let mut x = vec![0.0; n];
    for i in 0..n {
        x[i] = y[i] - factor * q[i];
    }

    x
}    