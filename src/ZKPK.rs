/* This file implements the ZKPK presented in https://eprint.iacr.org/2018/560 (pg.12-14)*/

use core::f64::consts::{E};
use crate::linearalgebra::{Vec1d,Vec2d,normal_vec2d,gauss_vec2d,Modq,Flatten,Norm,Transpose};
use crate::linearalgebraf::{power_iteration};
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
    pub n:usize,
    pub rho:f64,
}
#[derive(Debug,Clone)]
pub struct Proverst{
    /// *`A`: Fixed randomly choosen matrix
    /// *`S`: Preimages we want to prove we have knowdlege about 
    /// *`T`: Product of A*S
    /// *`T`: Standard deviation for gaussian distributions.
    /// *`B`: Bound of the norms of the matrix Z.
    pub A:Vec2d,
    pub T:Vec2d,
    pub S:Vec2d,
    pub sigma:f64,
    pub B:f64
}

pub struct Verifierst{
    /// *`A`: Fixed randomly choosen matrix
    /// *`T`: Product of A*S
    /// *`B`: Bound of the norms of the matrix Z.
    pub A:Vec2d,
    pub T:Vec2d,
    pub B:f64,
}
/// ZKPK constructor 
impl ZKPK{
    pub fn new(lambda:i64,q:i64,r:usize,l:usize,v:usize,rho:f64)->ZKPK{
        ZKPK{
            lambda:lambda,
            q:q,
            r:r, 
            l:l,
            v:v,
            n:(lambda+2) as usize, 
            rho:rho,
        }
    }
}
/// Proverst constructor
///
/// Arguments:
///
/// *`z`: ZKPK struct
/// *`S`:Preimages we want to prove we have knowdlege about 

impl Proverst{
    pub fn new(z:&ZKPK,S:&Vec2d) -> Proverst {
        
        let A=normal_vec2d(z.r,z.v,z.q);
        let sigma=((12.0)/(z.rho as f64).ln())*power_iteration(&(&S.transpose()*S),5)*(z.l as f64*z.n as f64).sqrt();
        Proverst {
            A: A.clone(),
            S:S.clone(),
            T:(&A*S).mod_q(z.q),
            sigma:sigma,
            B:(2.0*z.v as f64).sqrt()*2.0*sigma
        }
    }
}
/// Verifierst constructor
///
/// Arguments:
///
/// *`A`: ZKPK struct
/// *`T`: A*S
/// *`B`: Bound for the vector norms.

impl Verifierst{
    pub fn new(A:Vec2d,T:Vec2d,B:f64) -> Verifierst {

        Verifierst {
            A:A,
            T:T,
            B:B,
        }
    }
}
/// Definition of the trait Prover. It implements two functions, the computation of both W matrix and Z matrix.

pub trait Prover<T,U,K> {

    fn compute_w(&self,z:&U)->(T,T);
    fn compute_z(&self,c:&T,y:&T,rho:f64)->(T,K);

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
/// (Vec2d,bool). The first element of the tuple contains the computation of the matrix Z 
///while the second one is the result of the abort step.

impl Prover<Vec2d,ZKPK,bool> for Proverst{
    fn compute_w(&self,z:&ZKPK)->(Vec2d,Vec2d){    
        let Y=gauss_vec2d(0.0,self.sigma,6,120,z.v,z.n);
        ((&self.A*&Y),Y)
    }
    fn compute_z(&self,c:&Vec2d,Y:&Vec2d,rho:f64)->(Vec2d,bool){
        let Z=&(&self.S*c)+Y;
        let result=rejection_sampling(Z.clone(),&self.S*c,self.sigma,rho);
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
        let norm_ver=verify_norms(Z,self.B);
        let AZ=(&self.A*Z).mod_q(z.q);
        let TCW=(&(&self.T*C)+W).mod_q(z.q);
        (&AZ-&TCW).flatten().0.iter().sum::<i64>()==0 && norm_ver 
    }
}
/// Function to verify if the norms of each vector in the matrix Z are lower than B.
///
/// Arguments:
///
///* `Z`: Z matrix
/// * `bound`: B=2*sqrt(2*v)*sigma
///
///
/// Return:
///
///  Boolean: True if the previous condtion is fulfilled. False otherwise.
fn verify_norms(Z:&Vec2d,bound:f64)->bool{
    
    let ZT=Z.transpose();
    let norm_vec=(0..ZT.N).map(|i| Vec1d(ZT.matrix[i].clone()).norm()).collect::<Vec<f64>>();
    
    norm_vec.into_iter().find(| &x| x>=bound)==None
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
pub fn rejection_sampling(Z:Vec2d,B:Vec2d,sigma:f64,rho:f64)->bool{
    
    let flatten_z= Z.flatten();
    let flatten_b= B.flatten();
    
    let mut rng = rand::thread_rng();
    let u:f64=rng.gen::<f64>();
    
    let zbdot=(&flatten_z*&flatten_b) as f64;
    let bnorm=(&flatten_b*&flatten_b) as f64;
    
    let gauss_value=(1.0/rho)*E.powf((-2.0*zbdot+bnorm)/(2.0*((sigma).powf(2.0))));
    
    if u>gauss_value{
        return false;
    }
    
    else{
        return true;
    }
}