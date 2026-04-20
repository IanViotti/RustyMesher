#![allow(non_snake_case)]

use serde::Deserialize;

/// Main configuration structure that holds all parameters required for the 
/// grid generation process. This struct is populated by reading the `job_config.toml` file.
#[derive(Clone, Debug, Deserialize)]
pub struct Config {
    // Basic identifiers
    /// Name of the mesh/job, used for creating output directories and files.
    pub meshname: String,
    
    /// The physical shape of the inner boundary (e.g., Biconvex or NACA).
    pub airfoil_type: AirfoilType,  
    
    /// Outer boundary radius or distance. Defines how far the computational 
    /// domain extends from the airfoil surface.
    pub r_max: f64,
    
    // Grid topology parameters
    /// Number of points in the longitudinal direction (xi / ξ coordinate).
    /// Typically represents the wrap-around direction in an O-grid.
    pub longitudinal_points: usize,
    
    /// Number of points in the normal direction (eta / η coordinate).
    /// Represents the radial lines extending from the airfoil to the outer boundary.
    pub normal_points: usize,
    
    /// Geometric stretching factor used to cluster grid points closer to the 
    /// airfoil surface (viscous/boundary layer resolution) in the initial algebraic/parabolic grid.
    pub stretching_factor: f64, 
    
    // Geometry definition parameters
    /// Geometry parameter 'q', often used to define the specific shape or thickness 
    /// distribution of the biconvex airfoil.
    pub q: f64,
    
    /// Maximum thickness-to-chord ratio (t/c) of the airfoil.
    pub t: f64,
    
    // Numerical solver parameters
    /// Relaxation factor used in iterative schemes (like SOR - Successive Over-Relaxation).
    /// If omega = 1.0, it acts as standard Gauss-Seidel.
    pub omega: f64,
    
    /// Lower bound for the acceleration parameter sequence used in the 
    /// AF1 (Approximate Factorization) algorithm to speed up convergence.
    pub alpha_l: f64,
    
    /// Upper bound for the acceleration parameter sequence in the AF1 algorithm.
    pub alpha_h: f64,
    
    /// Selects the underlying partial differential equation (PDE) for the 
    /// elliptic grid generator (Laplace vs. Poisson with grid control).
    pub elliptic_equation: EllipticEquation, 
    
    /// Maximum number of iterations allowed for the elliptic solver before it stops.
    pub max_iter: usize,
    
    /// Convergence tolerance. The solver stops when the maximum residual 
    /// error falls below this value.
    pub conv_criterion: f64,
}

/// Defines the mathematical formulation used to generate the internal grid points.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
pub enum EllipticEquation {
    /// Solves the homogeneous Laplace equation (∇²ξ = 0, ∇²η = 0).
    /// Tends to produce smooth grids but lacks direct control over spacing near the boundaries.
    Laplace,
    
    /// Solves the inhomogeneous Poisson equation (∇²ξ = P, ∇²η = Q).
    /// Introduces source terms (P, Q) to enforce grid clustering and orthogonality near the boundaries.
    Poisson {
        /// Target normal grid spacing (ds) at the airfoil surface (inner boundary wall).
        ds_wall: f64, 
        
        /// Decay factor 'a' that dictates how quickly the longitudinal grid control 
        /// fades out as it moves away from the source.
        a_decay: f64, 
        
        /// Decay factor 'b' that dictates how quickly the normal grid control 
        /// (attraction to the wall) fades out towards the outer boundary.
        b_decay: f64, 
    }, 
}

/// Available airfoil geometries for the inner boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
pub enum AirfoilType {
    /// A symmetrical biconvex airfoil formed by circular arcs.
    Biconvex,
    
    /// A standard 4-digit symmetrical NACA airfoil (e.g., NACA 0012).
    NACA00XX, 
}