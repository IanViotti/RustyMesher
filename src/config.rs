#![allow(non_snake_case)]

use crate::geometry::AirfoilType;

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
    pub n_max: usize,
    pub conv_criterion: f64,
}