use core::f64::consts::{E};
use crate::linearalgebra::{Vec1d,Vec2d,normal_vec2d,gauss_vec2d,Modq,Flatten,Norm,Transpose};
use crate::linearalgebraf::{power_iteration};
use rand::Rng;
pub struct ZKPK{
    pub lambda:i64,
    pub q:i64,
    pub r:usize,
    pub l:usize,
    pub v:usize,
    pub n:usize,
    pub rho:i64,
}
#[derive(Debug,Clone)]
pub struct Proverst{
    pub A:Vec2d,
    pub T:Vec2d,
    pub S:Vec2d,
    pub sigma:f64,
    pub B:f64
}

pub struct Verifierst{
    pub A:Vec2d,
    pub T:Vec2d,
    pub B:f64,
}
impl ZKPK{
    pub fn new(lambda:i64,q:i64,r:usize,l:usize,v:usize,rho:i64)->ZKPK{
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
impl Proverst{
    pub fn new(z:&ZKPK,S:&Vec2d) -> Proverst {
        let A=normal_vec2d(z.r,z.v,z.q);
        let sigma=((12.0)/(z.rho as f64).ln())*power_iteration(S.clone(),5)*(z.l as f64*z.n as f64).sqrt();
        println!("sigma:{}",sigma);
        Proverst {
            A: A.clone(),
            S:S.clone(),
            T:(&A*S).mod_q(z.q),
            sigma:sigma,
            B:(2.0*z.v as f64).sqrt()*sigma
        }
    }
}
impl Verifierst{
    pub fn new(z:&ZKPK,A:Vec2d,T:Vec2d,B:f64) -> Verifierst {
        Verifierst {
            A:A,
            T:T,
            B:B,
        }
    }
}
pub trait Prover<T,U,K> {
    fn compute_w(&self,z:&U)->(T,T);
    fn compute_z(&self,c:&T,y:&T,z:&U,rho:f64)->(T,K);
}
pub trait Verifier<T,U> {
    fn compute_c(&self,z:&ZKPK)->T;
    fn verify(&self,Z:&T,C:&T,W:&T,z:&U)->bool;
}

impl Prover<Vec2d,ZKPK,bool> for Proverst{
    fn compute_w(&self,z:&ZKPK)->(Vec2d,Vec2d){
        let Y=gauss_vec2d(0.0,self.sigma,14,120,z.v,z.n);
        ((&self.A*&Y).mod_q(z.q),Y)
    }
    fn compute_z(&self,c:&Vec2d,Y:&Vec2d,z:&ZKPK,rho:f64)->(Vec2d,bool){
        let Z=(&(&self.S*c)+Y).mod_q(z.q);
        let result=rejection_sampling(Z.clone(),&self.S*c,self.sigma,rho);
        (Z,result)
    }
}
impl Verifier<Vec2d,ZKPK> for Verifierst{
    fn compute_c(&self,z:&ZKPK)->Vec2d{
        challenge_mat(z.l,z.n)
    }
    fn verify(&self,Z:&Vec2d,C:&Vec2d,W:&Vec2d,z:&ZKPK)->bool{
        let norm_ver=verify_norms(Z,self.B);
        let AZ=(&self.A*Z).mod_q(z.q);
        let TCW=(&(&self.T*C)+W).mod_q(z.q);
        println!("AZ:{:?}",AZ);
        println!("TCW:{:?}",TCW);
        (&AZ-&TCW).flatten().0.iter().sum::<i64>()==0 && norm_ver 
    }
}

fn verify_norms(Z:&Vec2d,bound:f64)->bool{
    let ZT=Z.transpose();
    let norm_vec=(0..ZT.N).map(|i| Vec1d(ZT.matrix[i].clone()).norm()).collect::<Vec<f64>>();
    norm_vec.into_iter().find(| &x| x>=bound)==None
}
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
fn rejection_sampling(Z:Vec2d,B:Vec2d,sigma:f64,rho:f64)->bool{
    let flatten_z= Z.flatten();
    let flatten_b= B.flatten();
    let mut rng = rand::thread_rng();
    let u:f64=rng.gen::<f64>();
    let zbdot=(&flatten_z*&flatten_b) as f64;
    let bnorm=(&flatten_b*&flatten_b) as f64;
    let gauss_value=(1.0/rho)*E.powf((-2.0*zbdot+bnorm)/2.0*((sigma).powf(2.0)));
    println!("zdbot:{}",zbdot);
    println!("bnorm:{}",bnorm);
    println!("u:{}",u);
    println!("gauss:{}",gauss_value);
    if u>gauss_value{
        return false;
    }
    else{
        return true;
    }
}