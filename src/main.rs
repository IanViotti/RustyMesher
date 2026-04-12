#![allow(non_snake_case)]

mod mesher_core;
mod mesher_utils;
mod geometry;
mod config;
mod parabolic_mesher;

fn main() {
    let config = config::Config {
        meshname: "biconvex_mesh".to_string(), // Mesh name for output organization
        r_max: 6.5,                   // Circunference radius of the outter mesh boundary
        longitudinal_points: 93,   // Number of points along the airfoil chord
        normal_points: 5,          // Number of points in the normal direction from the airfoil surface
        t: 0.10,                    // Airfoil thickness parameter
        stretching_factor: 1.2,     // Stretching factor for geometric progression (1.0 for uniform mesh)
        P: 1.0,                     // Fator de inclinação da malha
        Q: 1.1,                     // Fator de amortecimento
        n_max: 5000,                // Maximum number of solver iterations
        conv_criterion: 0.0,        // Convergence criterion for residual
    };

    println!("\n
            -----------------------------\n
            Starting mesh generation job: {}\n
            -----------------------------\n", config.meshname);

    // Create necessary directories for output
    mesher_utils::init_mesh_directory(&config);

    // Initialize griid
    let mut grid = mesher_utils::create_grid(config.longitudinal_points, config.normal_points);

    // Insert geometry into grid
    geometry::insert_geoemtry(&config, &mut grid);

    // create mesh
    mesher_core::create_mesh(&config, &mut grid);

    match mesher_utils::save_grid_to_csv(&grid, "job_files/biconvex_mesh/malha_gerada.csv") {
        Ok(_) => println!("Arquivo salvo com sucesso!\n"),
        Err(e) => eprintln!("Erro ao salvar o arquivo: {}", e),
    }

}
