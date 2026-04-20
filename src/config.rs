#![allow(non_snake_case)]

use crate::geometry::AirfoilType;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EllipticEquation {
    Laplace,
    Poisson {
        ds_wall: f64, // Espaçamento na parede
        a_decay: f64, // Decaimento longitudinal
        b_decay: f64, // Decaimento normal
    }, 
}


#[derive(Clone, Debug)]
pub struct Config {
    pub meshname: String,
    pub airfoil_type: AirfoilType,  
    pub r_max: f64,
    pub longitudinal_points: usize,
    pub normal_points: usize,
    pub stretching_factor: f64, // New parameter for geometric stretching
    pub q: f64,
    pub t: f64,
    pub omega: f64,
    pub alpha_l: f64,
    pub alpha_h: f64,
    pub elliptic_equation: EllipticEquation, // New parameter to select the elliptic equation type
    pub max_iter: usize,
    pub conv_criterion: f64,
}