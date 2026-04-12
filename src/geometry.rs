use crate::config;
use ndarray::Array2;
use std::f64::consts::PI;
use crate::mesher_utils::Point;


pub fn insert_geoemtry(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Inserting geometry...");

    // Insert the airfoil geometry into the grid
    insert_airfoil_geom(config, grid);

    // Insert the outer boundary geometry into the grid
    insert_outter_boundary(config, grid);

}


fn insert_airfoil_geom(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Inserting biconvex airfoil geometry...");

    let n_longitudinal = config.longitudinal_points;
    let t = config.t;
    
    // Supondo que você adicione 'stretching_factor' no seu config.rs
    // Um valor típico seria algo como 1.1 ou 1.05. Se for 1.0, a malha é uniforme.
    let r = config.stretching_factor; 

    // O número de pontos na superfície (de x=1 até x=0 e voltando para x=1)
    // Para facilitar, vamos pré-calcular a distribuição de x de 0.0 a 1.0 usando o stretching.
    // A quantidade de pontos entre o Bordo de Ataque (x=0) e o Bordo de Fuga (x=1) em uma face é:
    let n_face = (n_longitudinal - 1) / 2 + 1;
    let mut x_dist = vec![0.0; n_face];

    // Calculando a distribuição de x apenas para uma face (x=0 até x=1)
    // Agrupando no Bordo de Ataque (x=0) e Bordo de Fuga (x=1).
    // Faremos isso preenchendo do BA (0.0) até o meio (0.5), e espelhando para o resto.
    let n_half = (n_face - 1) / 2; // Número de intervalos até o meio da corda
    
    if r == 1.0 {
        // Malha uniforme
        for i in 0..n_face {
            x_dist[i] = i as f64 / (n_face - 1) as f64;
        }
    } else {
        // Distribuição com Progressão Geométrica
        // A soma da PG é S = a1 * (r^n - 1) / (r - 1). Queremos que S seja 0.5 (metade da corda)
        let a1 = 0.5 * (r - 1.0) / (r.powi(n_half as i32) - 1.0);
        
        x_dist[0] = 0.0;
        let mut current_dx = a1;
        
        for k in 1..=n_half {
            x_dist[k] = x_dist[k - 1] + current_dx;
            current_dx *= r; // Aumenta o dx pelo fator de estiramento
        }
        
        // Espelha a distribuição para a segunda metade (de x=0.5 até x=1.0)
        for k in 1..=n_half {
            let idx_mirror = n_half + k;
            x_dist[idx_mirror] = 1.0 - x_dist[n_half - k];
        }
    }

    // Agora, mapeamos essa distribuição de x_dist para o contorno do perfil na malha O.
    // i=0 começa no Bordo de Fuga (x=1.0), vai pelo intradorso até o Bordo de Ataque (x=0.0) em i = n_face - 1
    // e volta pelo extradorso até o Bordo de Fuga (x=1.0) em i = n_longitudinal - 1.
    for i in 0..n_longitudinal {
        let x = if i < n_face {
            // Intradorso (caminhando de x=1 para x=0)
            x_dist[n_face - 1 - i]
        } else {
            // Extradorso (caminhando de x=0 para x=1)
            x_dist[i - n_face + 1]
        };

        // Espessura teórica do perfil
        let thickness = 2.0 * t * x * (1.0 - x);

        // Define a coordenada y: intradorso (y negativo) ou extradorso (y positivo)
        let y = if i < n_face {
            -thickness 
        } else {
            thickness  
        };

        // Instanciamos o struct Point e salvamos na matriz computacional
        grid[[i, 0]] = Point { x, y };
    }
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