use crate::config;
use ndarray::Array2;
use crate::parabolic_mesher;
use crate::elliptic_mesher;
use crate::mesher_utils::Point;

pub fn create_mesh(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Initializing mesh generation process...");

    // Parabolic mesh generation 
    parabolic_mesher::generate_parabolic_mesh(config, grid);

    // Elliptic mesh generation 
    elliptic_mesher::elliptic_smoother(config, grid);

}