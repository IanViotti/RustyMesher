use crate::config;
use ndarray::Array2;
use crate::parabolic_mesher;
use crate::elliptic_mesher;
use crate::mesher_utils::Point;

/// Orchestrates the internal grid generation process.
/// In CFD (Computational Fluid Dynamics) mesh generation, it is common practice 
/// to first generate an initial internal point distribution using a fast, non-iterative 
/// method, and then refine it using a more robust, iterative elliptic solver.
pub fn create_mesh(config: &config::Config, grid: &mut Array2<Point>) {
    println!("Initializing mesh generation process...");

    // Step 1: Parabolic mesh generation (Marching Method)
    // Computes the initial guess for the internal grid coordinates. 
    // This method (often based on Nakamura's marching approach) propagates 
    // outward from the inner boundary (airfoil surface) towards the outer field.
    // It is computationally cheap and provides a very structured starting grid, 
    // which significantly reduces the number of iterations required by the elliptic smoother.
    parabolic_mesher::generate_parabolic_mesh(config, grid);

    // Step 2: Elliptic mesh generation (PDE Smoother)
    // Solves the specified Elliptic Partial Differential Equation (Laplace or Poisson) 
    // using an iterative numerical scheme (like ADI or AF1).
    // It takes the parabolic grid as its initial condition and spatially smooths it.
    // If the Poisson equation is selected, this step also introduces source terms (P, Q)
    // to enforce grid clustering and strict orthogonality at the solid boundaries 
    // (commonly utilizing Sorenson/Steger grid control techniques).
    elliptic_mesher::elliptic_smoother(config, grid);

}