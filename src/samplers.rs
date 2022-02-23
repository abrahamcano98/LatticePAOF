/* This file contains the sampling of discrete gaussian distributions as presented in https://ieeexplore.ieee.org/document/8314133 */
#![allow(non_snake_case)]

use rand::{distributions::Uniform, Rng};
use std::collections::HashMap;
use core::f64::consts::{E,PI};

/// Function to sample a uniformely distributed vector over Zq
/// Arguments:
///
/// * `size`: Lenght of the vector
/// * `max_value`: Q that generates the field
///
///
/// Return:
///
///  Vec<i64> with the normal vector sampled

pub fn normal_vec(size:usize, max_value:i64)->Vec<i64>{
    let mut rng = rand::thread_rng();
    let range = Uniform::new(0, max_value);
    (0..size).map(|_| rng.sample(&range)).collect::<Vec<i64>>()
}

/// Function to sample a uniformely distributed matrix over Zq
/// Arguments:
///
/// * `size`: Lenght of the vector
/// * `max_value`: Q that generates the field
///
///
/// Return:
///
///  Vec<Vec<i64>> with the normal matrix sampled.
pub fn normal_mat(height:usize, width:usize, max_value:i64)->Vec<Vec<i64>>{
    let mut rng = rand::thread_rng();
    let range = Uniform::new(0, max_value);
    (0..height)
        .map(|_| {
            (0..width)
                .map(|_| rng.sample(&range)).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
}
/// Function to get the integer and decimal part given a number.
/// Arguments:
///
/// * `x`: Number to be split
/// * `precision`: Number of decimals to be obtained.
///
///
/// Return:
///
///  (Y,Z) where Y represents the integer part and Z the decimal part.
fn modf(x:f64,precision:f64)->(i64,i64){
    let y= x as i64;
    let mut z=x-(y as f64);
    z=z*precision;
    z=(z*precision).round()/precision;
    (y,z as i64)
}
/// Function to obtain the binary representation of a float number
/// Arguments:
///
/// * `num`: Number which we wish to convert.
/// * `precision`: Lenght of the binary decimal. The larger, the more accurate.
///
/// Return:
///
///  Vec<i64>, whose coefficients represent the binary representation.

fn float_to_binary(mut num:f64,precision:usize)->Vec<i64>{
    let mut binary_num=vec![0i64;precision];
    for i in 0..precision{
        num=num*2.0;
        let (integer,_)=modf(num,1.0);
        if integer==1{
            binary_num[i]=1;
            num=num-1.0;
        } else{
            binary_num[i]=0;
        }
    }
    binary_num
}
/// Function to obtain the gaussian probability of a number
/// Arguments:
///
/// * `mu`: Mean of the gaussian.
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// Return:
///
///  Hashmap<i64,i64>, Vec<f64> The hashmap has the equivalence between the values computed and the index in the vector. The vector contains the probabilities itself.
fn compute_gaussian_probability(mu:f64,sigma:f64,t:i64)->(HashMap<i64,i64>,Vec<f64>){
   
    let max_range=t*(sigma as i64);
    let gauss_map:HashMap<i64,i64>=(0..2*max_range).zip(-max_range..max_range).collect();
    let mut prob_vec=vec![0f64;(2*max_range).try_into().unwrap()];


    for i in 0..2*max_range{
        
        prob_vec[i as usize]=(1.0/((sigma)*((2.0*PI).sqrt())))*(E.powf((-(gauss_map[&i] as f64-mu).powf(2.0))/(2.0*sigma.powf(2.0))))
    }
    (gauss_map,prob_vec)
}
/// Function to sample a matrix, whose rows contain the binary expresions of the probabilities of the integers between (-t*sigma..t*sigma)
///
/// Arguments:
///
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// * `precision`: Gaussian probabilities decimals. The greater the number the more approximate the gaussian.
/// * `prob_vec`: Probability vector
///
///
/// Return:
///
///  Vec<Vec<i64>> with the matrix sampled

fn binary_prob_mat(sigma:f64,t:i64,precision:usize,prob_vec:Vec<f64>)->Vec<Vec<i64>>{
    let max_range=t*(sigma as i64);
    let mut knuth_mat=Vec::new();
    //println!("Prob_vec:{:?}",prob_vec);
    for i in 0..2*max_range{
        knuth_mat.push(float_to_binary(prob_vec[i as usize],precision))
    }
    knuth_mat
}
/// Function to get a single value accoding to a gaussian distribution with relation to the knuth-yao algorithm.
///
/// Arguments:
///
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// * `gauss_map`: Hashmap with the equivalence between indexes and integers between (-t*sigma..t*sigma)
/// * `knuth-mat`: Binary probability matrix
///
///
/// Return:
///
///  i64: Value sampled.
fn knuth_yao(sigma:f64,t:i64,gauss_map:HashMap<i64,i64>,knuth_mat:Vec<Vec<i64>>)->i64{
    let max_range=t*(sigma as i64);
    let mut d=0;
    let mut col=0;
    let mut rng = rand::thread_rng();
    loop{        
        let r: i64 = rng.gen::<f64>().round() as i64;
        d=2*d+r;
        for row in (0..2*max_range).rev(){
            d=d-knuth_mat[row as usize][col as usize];
            if d==-1{
                return gauss_map[&row];
            }
        }
        col=col+1;
    }
}
/// Function to get a guassian vector accoding to a gaussian distribution by means of the knuth-yao algorithm.
///
/// Arguments:
///
///* `mu`: Mean of the gaussian distribution.
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// * `precision`: Number of decimals for the gaussian probabilities.
/// * `length`: Length of the vector.
///
///
/// Return:
///
///  Vec<i64> with the vector sampled
pub fn gauss_vec(mu:f64,sigma:f64,t:i64, precision:usize, length:usize)->Vec<i64>{
    let (gauss_map,prob_vec)=compute_gaussian_probability(mu,sigma,t);
    let knuth_mat=binary_prob_mat(sigma,t,precision,prob_vec.clone()); 
    let gauss_vec:Vec<i64>=(0..length).map(|_| knuth_yao(sigma,t,gauss_map.clone(),knuth_mat.clone())).collect();
    return gauss_vec;
}
/// Function to get a guassian matrix accoding to a gaussian distribution by means of the knuth-yao algorithm.
///
/// Arguments:
///
///* `mu`: Mean of the gaussian distribution.
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// * `precision`: Number of decimals for the gaussian probabilities.
/// * `length`: Rows of the matrix.
/// *`width`: Columns of the matrix
///
///
/// Return:
///
///  Vec<i64> with the matrix sampled
pub fn gauss_mat(mu:f64,sigma:f64,t:i64, precision:usize, length:usize,width:usize)->Vec<Vec<i64>>{
    let (gauss_map,prob_vec)=compute_gaussian_probability(mu,sigma,t);
    let knuth_mat=binary_prob_mat(sigma,t,precision,prob_vec.clone()); 
    let gauss_mat:Vec<Vec<i64>>=(0..length)
    .map(|_| {
        (0..width)
            .map(|_| knuth_yao(sigma,t,gauss_map.clone(),knuth_mat.clone())).collect::<Vec<i64>>()
    })
    .collect::<Vec<Vec<i64>>>();
    return gauss_mat;
}