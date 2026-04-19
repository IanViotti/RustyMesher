#![allow(non_snake_case)]

use std::time::Instant;
mod mesher_core;
mod mesher_utils;
mod geometry;
mod config;
mod parabolic_mesher;
mod elliptic_mesher;
use crate::geometry::AirfoilType::Biconvex;
use crate::geometry::AirfoilType::NACA00XX;

fn main() {

    let start = Instant::now();

    let config = config::Config {
        meshname: "naca0012_mesh".to_string(), // Mesh name for output organization
        airfoil_type: NACA00XX,
        r_max: 6.5,                   // Circunference radius of the outter mesh boundary
        longitudinal_points: 93,   // Number of points along the airfoil chord
        normal_points: 15,          // Number of points in the normal direction from the airfoil surface
        t: 0.12,                    // Airfoil thickness parameter
        stretching_factor: 1.2,     // Stretching factor for airfoil surface point distribution (1.0 for uniform mesh)
        q: 1.15,                     // Stretching factor for the parabolic mesher (1.0 for no stretching)
        omega: 1.5,                 // Relaxation factor
        alpha_l: 0.5,               // Minimum alpha for multi-frequency acceleration
        alpha_h: 10.0,              // Maximum alpha for multi-frequency acceleration
        max_iter: 5000,                // Maximum number of solver iterations
        conv_criterion: 1e-6,        // Convergence criterion for residual
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

    let mesh_filepath = format!("job_files/{}/mesh.csv", config.meshname);
    match mesher_utils::save_grid_to_csv(&grid, &mesh_filepath) {
        Ok(_) => println!("File saved as .csv..."),
        Err(e) => eprintln!("Erro ao salvar o arquivo: {}", e),
    }

    // Export to VTK for visualization
    let mesh_filepath = format!("job_files/{}/mesh.vtk", config.meshname);
    _ = mesher_utils::export_vtk_structured_grid(&grid, &mesh_filepath);

    
    println!("Job completed!");
    
    let duration = start.elapsed();
    println!("Tempo total da geração: {:?} \n", duration);

}
