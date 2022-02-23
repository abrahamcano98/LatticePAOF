/* This file implements the ZKPK presented in https://eprint.iacr.org/2018/560 (pg.12-14)*/

#![allow(non_snake_case)]

use core::f64::consts::{E};
use crate::linearalgebra::{Vec2d,normal_vec2d,gauss_vec2d,Modq,Flatten};
use rand::Rng;

/// Struct ZKPK which has the parameters used by the scheme.
pub struct ZKPK{
    /// *`lambda`: represents the security parameter
    /// *`q`: The generator of the field
    /// *`r`: Size of the commitments
    /// *`l`: Columns of the secret matrix
    /// *`v`: Rows of the secret matrix
    /// *`sigma`: Standard deviation of the error distribution
    /// *`n`: Columns of the challenge matrix
    pub lambda:i64,
    pub q:i64,
    pub r:usize,
    pub l:usize,
    pub v:usize,
    pub sigma:f64,
    pub n:usize,
}
/// Struct with the prover fields.
#[derive(Debug,Clone)]
pub struct Proverst{
    /// *`A`: Fixed randomly choosen matrix
    /// *`S`: Preimages we want to prove we have knowdlege about 
    /// *`T`: Product of A*S
    pub A:Vec2d,
    pub T:Vec2d,
    pub S:Vec2d,
}
/// Struct with the prover fields.
pub struct Verifierst{
    /// *`A`: Fixed randomly choosen matrix
    /// *`T`: Product of A*S
    pub A:Vec2d,
    pub T:Vec2d
}
/// Proverst constructor
///
/// Parameters:
///
/// *`z`: ZKPK struct
/// *`S`:Preimages we want to prove we have knowdlege about 
impl Proverst{
    pub fn new(z:&ZKPK,S:&Vec2d) -> Proverst {
        let A=normal_vec2d(z.r,z.v,z.q);
        Proverst {
            A: A.clone(),
            S:S.clone(),
            T:(&A*S).mod_q(z.q),
        }
    }
}
/// Verifierst constructor
impl Verifierst{
    pub fn new(p:&Proverst) -> Verifierst {
        Verifierst {
            A:p.clone().A,
            T:p.clone().T,
        }
    }
}
/// Definition of the trait Prover. It implements two functions, the computation of both W matrix and Z matrix.

pub trait Prover<T,U,K> {
    fn compute_w(&self,z:&U)->(T,T);
    fn compute_z(&self,c:&T,y:&T,z:&U,rho:f64)->(T,K);
}
/// Definition of the trait Verifier. It implements the sample of the challenge matrix C as well as the the verification of the proof.
pub trait Verifier<T,U> {
    fn compute_c(&self,z:&ZKPK)->T;
    fn verify(&self,Z:&T,C:&T,W:&T,z:&U)->bool;
}
/// Implementation of Prover for Provest.
///
/// Functions: compute_w, compute_z
///
/// fn compute_w:
///
/// Parameters:
///
///*z: ZKPK struct
///
/// return:
///
/// (Vec2d,Vec2d) where the first element of the tuple is the matrix W and the second one the matrix Y 
///
///fn compute_z:
///
/// Parameters:
///
///*c: Challenge matrix.
///*y: Gaussian matrix sampled by the prover. 
///*z: ZKPK struct
///*rho: Parameter for the artificial abort. It should abourt 1/rho times.
///
/// return
///
/// (Vec2d,bool). The first element of the tuple contains the computation of the matrix Z while the second one the result of the abort step.
impl Prover<Vec2d,ZKPK,bool> for Proverst{
    fn compute_w(&self,z:&ZKPK)->(Vec2d,Vec2d){
        let Y=gauss_vec2d(0.0,z.sigma,14,120,z.v,z.n);
        ((&self.A*&Y).mod_q(z.q),Y)
    }
    fn compute_z(&self,c:&Vec2d,Y:&Vec2d,z:&ZKPK,rho:f64)->(Vec2d,bool){
        let Z=(&(&self.S*c)+Y).mod_q(z.q);
        let result=rejection_sampling(Z.clone(),&self.S*c,z.sigma,rho);
        (Z,result)
    }
}
/// Implementation of Verifier for Verifierst.
///
/// Functions: compute_c, verify
///
/// fn compute_w:
///
/// Parameters:
///
///*z: ZKPK struct
///
/// return:
///
/// Vec2d Challenge matrix
///
///fn verify:
///
/// Parameters:
///
///*c: Challenge matrix.
///*w: W matrix computed by the prover. 
///*Z: Matrix Z computed by the prover.
///*z: ZKPK struct.

impl Verifier<Vec2d,ZKPK> for Verifierst{
    fn compute_c(&self,z:&ZKPK)->Vec2d{
        challenge_mat(z.l,z.n)
    }
    fn verify(&self,Z:&Vec2d,C:&Vec2d,W:&Vec2d,z:&ZKPK)->bool{
        let AZ=(&self.A*Z).mod_q(z.q);
        let TCW=(&(&self.T*C)+W).mod_q(z.q);
        (&AZ-&TCW).flatten().0.iter().sum::<i64>()==0
    }
}
/// Function to sample a challenge matrix. It is a matrix whose coefficients are either 0 or 1.
///
/// Arguments:
///
///* `height`: Rows of the matrix
/// * `width`: Columns of the matrix
///
///
/// Return:
///
///  Vec2d Challenge matrix.
fn challenge_mat(height:usize, width:usize)->Vec2d{
    let mut rng = rand::thread_rng();
    Vec2d{
        N:height,
        M:width,
        matrix:(0..height)
        .map(|_| {
            (0..width)
                .map(|_| rng.gen::<f64>().round() as i64).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
    }
}
/// Function to perform the artificial abort. 
///
/// Arguments:
///
/// * `Z`: Z matrix computed by the prover
/// * `B`: S*C matrix
/// * `sigma`: Standard deviation of the error distribution.
/// * `rho`: Control parameter. It should abort 1/rho times.
/// Return:
///
///  Boolean. True if abort, false otherwise.
fn rejection_sampling(Z:Vec2d,B:Vec2d,sigma:f64,rho:f64)->bool{
    let flatten_z= Z.flatten();
    let flatten_b= B.flatten();
    let mut rng = rand::thread_rng();
    let u:f64=rng.gen::<f64>();
    let a=-2.0*(&flatten_z*&flatten_b) as f64;
    let b=(&flatten_b*&flatten_b) as f64;
    let c=2.0*((sigma as f64).powf(2.0));
    let gauss_value=(1.0/rho)*E.powf((a+b)/c);
    if u>gauss_value{
        return false;
    }
    else{
        return true;
    }
}