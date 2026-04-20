use std::time::Instant;
use std::fs;
use toml;

// Import internal modules responsible for different aspects of the grid generation process.
mod mesher_core;      // Core logic for coordinating the mesh generation flow
mod mesher_utils;     // Utility functions (I/O operations, grid initialization, directory management)
mod geometry;         // Defines the physical boundaries (e.g., the biconvex airfoil surface)
mod config;           // Handles parsing and storing the TOML configuration parameters
mod parabolic_mesher; // Implements parabolic grid generation techniques (marching methods)
mod elliptic_mesher;  // Implements elliptic grid generation (Poisson/Laplace equations via ADI/AF1)

fn main() {

    // Start a timer to track the total execution time of the mesh generation
    let start = Instant::now();

    // 1. Read the TOML configuration file
    // This file contains user-defined parameters such as grid dimensions (longitudinal/normal points),
    // clustering/stretching factors, and the chosen numerical scheme (elliptic or parabolic).
    let config_content = fs::read_to_string("job_config.toml")
        .expect("It was not possible to read the job_config.toml file. Please make sure it exists and is in the correct path.\n");

    // 2. Deserialize using the 'toml' crate
    // Maps the raw string into a strongly-typed `config::Config` struct.
    let config: config::Config = toml::from_str(&config_content)
        .expect("Error parsing the job_config.toml file. Please check the format and values.\n");

    // Display job startup information
    println!("\n
            -----------------------------\n
            Starting mesh generation job: {}\n
            -----------------------------\n", config.meshname);

    println!("Configuration loaded: \n{:#?}\n", config);

    // Create necessary directories for output
    // Ensures that the destination folders for CSV and VTK files exist (e.g., job_files/meshname/)
    mesher_utils::init_mesh_directory(&config);

    // Initialize grid
    // Allocates the initial data structure for the computational domain (xi, eta)
    // based on the number of points in the longitudinal and normal directions.
    let mut grid = mesher_utils::create_grid(config.longitudinal_points, config.normal_points);

    // Insert geometry into grid
    // Applies the inner boundary conditions, which usually defines the 
    // physical shape of the biconvex airfoil at the lowest 'eta' level.
    geometry::insert_geoemtry(&config, &mut grid);

    // Create mesh
    // Dispatches the grid to the selected mathematical solver (elliptic or parabolic)
    // to compute the internal grid points coordinates in the physical space.
    mesher_core::create_mesh(&config, &mut grid);

    // Define path and save the final computed grid coordinates to a CSV file.
    // Useful for simple tabular data inspection and custom post-processing scripts.
    let mesh_filepath = format!("job_files/{}/mesh.csv", config.meshname);
    match mesher_utils::save_grid_to_csv(&grid, &mesh_filepath) {
        Ok(_) => println!("File saved as .csv..."),
        Err(e) => eprintln!("Erro ao salvar o arquivo: {}", e),
    }

    // Export to VTK for visualization
    // Saves the structured grid in a VTK format, which is the standard format for robust 
    // CFD visualization tools like ParaView or VisIt.
    let mesh_filepath = format!("job_files/{}/mesh.vtk", config.meshname);
    _ = mesher_utils::export_vtk_structured_grid(&grid, &mesh_filepath);

    // Print final completion message and elapsed execution time.
    println!("Job completed!");
    
    let duration = start.elapsed();
    println!("Total running time: {:?} \n", duration);

}