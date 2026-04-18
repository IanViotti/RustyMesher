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
    let r = config.stretching_factor; 

    // O índice do Bordo de Ataque será exatamente a metade do vetor.
    // IMPORTANTE: Para que a simetria seja perfeita, n_longitudinal DEVE ser ímpar (ex: 101).
    let mid = (n_longitudinal - 1) / 2;
    
    // Quantidade de pontos únicos que descrevem do Bordo de Ataque (0.0) ao Bordo de Fuga (1.0)
    let num_pts = mid + 1; 
    let mut x_dist = vec![0.0; num_pts];

    // 1. Gera a distribuição de x apenas de 0.0 a 1.0
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
        
        // Preenche agrupando de x=0.0 até x=0.5
        for k in 1..=n_half {
            x_dist[k] = x_dist[k - 1] + current_dx;
            current_dx *= r; 
        }
        
        // Espelha agrupando de x=0.5 até x=1.0
        for k in (n_half + 1)..num_pts {
            let mirror_idx = (num_pts - 1) - k;
            x_dist[k] = 1.0 - x_dist[mirror_idx];
        }
    }

    // 2. Mapeia a distribuição gerada para o contorno em anel da Malha O
    for i in 0..n_longitudinal {
        // Define a coordenada x descendo e subindo o array x_dist
        let x = if i <= mid {
            // Intradorso: começa em i=0 (pega x_dist[mid] = 1.0) 
            // Vai descendo até i=mid (pega x_dist[0] = 0.0)
            x_dist[mid - i]
        } else {
            // Extradorso: começa em i=mid+1 (pega x_dist[1])
            // Vai subindo até i=n_long-1 (pega x_dist[mid] = 1.0)
            x_dist[i - mid]
        };

        // Calcula a espessura teórica do perfil biconvexo
        let thickness = 2.0 * t * x * (1.0 - x);

        // Define a coordenada y
        let y = if i <= mid {
            -thickness // Intradorso (varredura de i=0 até i=mid)
        } else {
            thickness  // Extradorso (varredura de i=mid até o final)
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