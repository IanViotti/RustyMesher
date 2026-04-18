use crate::config;
use ndarray::Array2;
use std::f64::consts::PI;
use crate::mesher_utils::Point;

// Adicione este enum no seu arquivo config.rs ou no próprio módulo do mesher
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AirfoilType {
    Biconvex,
    NACA00XX, // Nomeado genericamente, pois a espessura 't' define se é 0012, 0015, etc.
}


pub fn insert_geoemtry(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Inserting geometry...");

    // Insert the airfoil geometry into the grid
    // Redireciona para a função correta baseada no enum configurado
    match config.airfoil_type {
        AirfoilType::Biconvex => insert_biconvex_geom(config, grid),
        AirfoilType::NACA00XX => insert_naca00xx_geom(config, grid),
    }

    // Insert the outer boundary geometry into the grid
    insert_outter_boundary(config, grid);

}

fn insert_biconvex_geom(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Generating Biconvex airfoil geometry...");

    let n_longitudinal = config.longitudinal_points;
    let t = config.t;
    let mid = (n_longitudinal - 1) / 2;
    
    // Chama a função auxiliar para pegar as coordenadas x
    let x_dist = generate_x_distribution(config);

    for i in 0..n_longitudinal {
        let x = if i <= mid { x_dist[mid - i] } else { x_dist[i - mid] };
        let x_safe = x.max(0.0);

        let thickness = 2.0 * t * x_safe * (1.0 - x_safe);

        let y = if i <= mid { -thickness } else { thickness };

        grid[[i, 0]] = Point { x, y };
    }
}

fn insert_naca00xx_geom(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Generating NACA 00xx (Left-Pointing: Nose at X=1, Tail at X=0)...");

    let n_longitudinal = config.longitudinal_points;
    let t = config.t;
    let mid = (n_longitudinal - 1) / 2;
    
    // Vetor com a distribuição de pontos matemáticos (de 0.0 até 1.0)
    let x_dist = generate_x_distribution(config);

    for i in 0..n_longitudinal {
        // x_math dita a posição ao longo da corda. 
        // 1.0 é a cauda, 0.0 é o nariz.
        let x_math = if i <= mid {
            x_dist[mid - i] // Descendo da cauda (1.0) para o nariz (0.0)
        } else {
            x_dist[i - mid] // Subindo do nariz (0.0) de volta para a cauda (1.0)
        };

        let x_safe = x_math.max(0.0);

        // Equação NACA 4 dígitos (com fechamento afiado em -0.1036 para o bordo de fuga)
        let thickness = 5.0 * t * (
            0.2969 * x_safe.sqrt() - 
            0.1260 * x_safe - 
            0.3516 * x_safe.powi(2) + 
            0.2843 * x_safe.powi(3) - 
            0.1036 * x_safe.powi(4)
        );

        // === A ÚNICA MUDANÇA PARA APONTAR PARA A ESQUERDA ===
        // Quando x_math for 1.0 (cauda), X físico será 0.0
        // Quando x_math for 0.0 (nariz), X físico será 1.0
        let x_coord = x_safe;

        // Mantém a varredura rigorosamente como você pediu:
        // i <= mid: Vai por baixo (Intradorso)
        // i > mid: Volta por cima (Extradorso)
        let y_coord = if i <= mid {
            -thickness 
        } else {
            thickness  
        };

        grid[[i, 0]] = Point { x: x_coord, y: y_coord };
    }
}

/// Gera a distribuição de pontos em X (0.0 a 1.0) com agrupamento (stretching)
fn generate_x_distribution(config: &config::Config) -> Vec<f64> {
    let n_longitudinal = config.longitudinal_points;
    let r = config.stretching_factor; 

    let mid = (n_longitudinal - 1) / 2;
    let num_pts = mid + 1; 
    let mut x_dist = vec![0.0; num_pts];

    if r <= 1.0 {
        // Malha uniforme
        for k in 0..num_pts {
            x_dist[k] = (k as f64) / ((num_pts - 1) as f64);
        }
    } else {
        // Distribuição com Progressão Geométrica
        let n_half = (num_pts - 1) / 2; 
        let a1 = 0.5 * (r - 1.0) / (r.powi(n_half as i32) - 1.0);
        
        x_dist[0] = 0.0;
        let mut current_dx = a1;
        
        for k in 1..=n_half {
            x_dist[k] = x_dist[k - 1] + current_dx;
            current_dx *= r; 
        }
        
        for k in (n_half + 1)..num_pts {
            let mirror_idx = (num_pts - 1) - k;
            x_dist[k] = 1.0 - x_dist[mirror_idx];
        }
    }

    x_dist
}

fn insert_outter_boundary(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Inserting outer geometry...");

    let n_longitudinal = config.longitudinal_points;
    
    // Supondo que você adicionou r_max (Raio do domínio externo) no config.
    // Se não, você pode fixar um valor aqui (ex: 20.0).
    let r_max = config.r_max; 

    // Obtém a dimensão normal (eta) diretamente do tamanho da matriz grid
    let (_, n_normal) = grid.dim();
    let j_max = n_normal - 1; // Índice da fronteira externa

    // Centro do círculo externo (meio da corda do perfil biconvexo)
    let x_center = 0.5;
    let y_center = 0.0;

    for i in 0..n_longitudinal {
        // Varredura angular de 0 a -2*PI.
        // Isso garante que o sentido acompanhe a fronteira interna:
        // i = 0 (TE, em baixo) -> ângulo 0
        // i = n_face (LE)      -> ângulo -PI
        // i = final (TE, cima) -> ângulo -2*PI
        let theta = -2.0 * PI * (i as f64) / ((n_longitudinal - 1) as f64);

        // Equação paramétrica do círculo
        let x = x_center + r_max * theta.cos();
        let y = y_center + r_max * theta.sin();

        // Instancia o struct Point e aloca na última linha da matriz (j_max)
        grid[[i, j_max]] = Point { x, y };
    }
}