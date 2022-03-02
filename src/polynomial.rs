/* This file contains the necessary traits to performs polynomial operations over the ring Zq[x]/<x^n+1>*/
#![allow(non_snake_case)]

use std::ops::{Mul,Add,Sub};
use crate::linearalgebra::{Vec1d,Vec2d,Modq as Modqla,Pad,circulant_mat,MattoVec,normal_vec1d,gauss_vec1d,Transpose};


///Polynomial is the struct that represents the polynomials over the finite field Zq[x]/<x^n+1>

#[derive(Debug,Clone)]
pub struct Polynomial{
    ///n: Length of the field
    ///q: q that generates the field
    ///coeffs: coeffs of a polynomial.

    pub n:i64,
    pub q:i64,
    pub coeffs:Vec1d,}
/// Implementation of the struct Polynomial vector
#[derive(Debug,Clone)]
pub struct PolynomialVec{
    ///n: Length of the field
    ///q: q that generates the field
    ///vec: Vector of polynomials.
    
    pub n:i64,
    pub q:i64,
    pub vec:Vec<Polynomial>}

/// Implementation of polynomial constructor.
///
/// Basic usage:
/// let p=Polynomial::new(128,4099,Vec1d(vec![1,0,0]))
impl Polynomial{
    pub fn new(n:i64, q:i64, coeffs:Vec1d) -> Polynomial {
        Polynomial {
            n: n,
            q: q,
            coeffs:coeffs.pad(0,n as usize,0)
        }
    }
}
/// Implementation of the PolynomialVec constructor.
///
/// Basic usage:
/// let mut vec=Vec::new();
/// let p1=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let p2=Polynomial::new(128,4099,Vec1d(vec![0,1,0]));
/// vec.push(p1);
/// vec.push(p2);
/// let p_vec=PolynomialVec::new(128,4099,vec);

impl PolynomialVec{
    pub fn new(n:i64, q:i64, vec:Vec<Polynomial>) -> PolynomialVec {
        PolynomialVec {
            n: n,
            q: q,
            vec:vec
        }
    }
}

/// Definition of trait Modq for polynomials.
pub trait Modq{
    fn mod_q(self)-> Self ;
}
/// Implementation of mod q for polynomial. It computes the mod q of each coefficient.
///
/// Basic usage:
///
/// let p=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let pq=p.mod_q();
impl Modq for Polynomial{
    fn mod_q(mut self)-> Self {
        self.coeffs=self.coeffs.mod_q(self.q);
        self
    }
}
/// Rust build-in trait to perform the addition of two Polynomial instances 
///
///
////// Arguments:
///
/// V:&Polynomial: Reference to a Polynomial instance.
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p1=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let p1=Polynomial::new(128,4099,Vec1d(vec![1,1,0,3,4,6]));
/// let p3=&p1+&p2

impl Add<&Polynomial> for &Polynomial
{
    type Output = Polynomial;

    fn add(self, v: &Polynomial) -> Polynomial  {
        Polynomial::new(self.n,self.q,&self.coeffs+ &(v.coeffs))
    }
}
/// Rust build-in trait to perform the substraction of two Polynomial instances 
///
///
////// Arguments:
///
/// V:&Polynomial: Reference to a Polynomial instance.
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p1=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let p1=Polynomial::new(128,4099,Vec1d(vec![1,1,0,3,4,6]));
/// let p3=&p1-&p2

impl Sub<&Polynomial> for &Polynomial
{
    type Output = Polynomial;

    fn sub(self, v: &Polynomial) -> Polynomial  {
        Polynomial::new(self.n,self.q,&self.coeffs- &v.coeffs)
    }
}
/// Rust build-in trait to perform the multiplication of a Polynomial by a scalar. 
///
///
////// Arguments:
///
/// V:&Polynomial: Reference to a Polynomial instance.
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let p5=5*&p1

impl Mul<&Polynomial> for i64{
    type Output=Polynomial;
    
    fn mul(self,v:&Polynomial)->Polynomial{
        Polynomial::new(v.n,v.q,self*&v.coeffs)
        
    }
}
/// Rust build-in trait to perform the multiplication of a Polynomial by a scalar. 
///
///
////// Arguments:
///
/// V:i64: Scalar
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let p5=&p1*5

impl Mul<i64> for &Polynomial{
    type Output=Polynomial;
    
    fn mul(self, v:i64)-> Polynomial{
        Polynomial::new(self.n,self.q,v*&self.coeffs)
    }
}
/// Rust build-in trait to perform the multiplication of two polynomails by means of the Toeplitz Matrix Vector Product.
///
///
////// Arguments:
///
/// V:&Polynomial: Reference to a Polynomial instance.
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p1=Polynomial_vec::new(128,4099,Vec1d(vec![1,0,0]));
/// let p2=Polynomial::new(128,4099,Vec1d(vec![1,1,0,3,4,6]));
/// let p3=&p1*&p2

impl Mul<&Polynomial> for &Polynomial{
    type Output=Polynomial;
    fn mul(self, v:&Polynomial)->Polynomial{
        Polynomial::new(self.n,self.q,(&v.coeffs*&circulant_mat(&self.coeffs)).mattovec())
    }
}
/// Rust build-in trait to perform the multiplication of a polynomial_vec 
///by a scalar (element-wise)
///
///
////// Arguments:
///
/// V:&PolynomialVec: Reference to a PolynomialVec instance.
///
///
/// Return:
///
/// PolynomialVec instance.
///
/// Basic usage:
///
/// let p_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p5=5&p_v;

impl Mul<&PolynomialVec> for i64{
    type Output=PolynomialVec;
    fn mul(self, v:&PolynomialVec)->PolynomialVec{
        let vec=v.vec.iter().map(|e| (self*e)).collect::<Vec<Polynomial>>();
        PolynomialVec{n:v.n,q:v.q,vec}
    }
}
/// Rust build-in trait to perform the multiplication of a polynomial_vec 
///by a scalar (element-wise)
///
///
////// Arguments:
///
/// V:&PolynomialVec: Reference to a PolynomialVec instance.
///
///
/// Return:
///
/// PolynomialVec instance.
///
/// Basic usage:
///
/// let p_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p5=&p_v*5;
impl Mul<i64> for &PolynomialVec{
    type Output=PolynomialVec;
    fn mul(self, v:i64)->PolynomialVec{
       let vec=self.vec.iter().map(|e| (e*v)).collect::<Vec<Polynomial>>();
       PolynomialVec{n:self.n,q:self.q,vec:vec}
    }
}
/// Rust build-in trait to perform the addition of two polynomial_vec 
///instances (element-wise)
///
///
////// Arguments:
///
/// V:&PolynomialVec: Reference to a PolynomialVec instance.
///
///
/// Return:
///
/// PolynomialVec instance.
///
/// Basic usage:
///
/// let p1_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p2_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p3=&p1+&p2;
impl Add<&PolynomialVec> for &PolynomialVec{
    type Output=PolynomialVec;
    fn add(self, v:&PolynomialVec)-> PolynomialVec{
        if self.vec.len()!=v.vec.len(){
            panic!("Mismatched vectors size.")
        }
        let vec=self.vec.iter().zip(v.vec.iter())
        .map(|(b, v)| Polynomial::new(self.n,self.q,&b.coeffs + &v.coeffs))
        .collect::<Vec<Polynomial>>();
        PolynomialVec{n:self.n,q:self.q,vec:vec}
    }
}
/// Rust build-in trait to perform the substraction of two polynomial_vec 
///instances (element-wise)
///
///
////// Arguments:
///
/// V:&PolynomialVec: Reference to a PolynomialVec instance.
///
///
/// Return:
///
/// PolynomialVec instance.
///
/// Basic usage:
///
/// let p1_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p2_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p3=&p1-&p2;
impl Sub<&PolynomialVec> for &PolynomialVec{
    type Output=PolynomialVec;
    fn sub(self, v:&PolynomialVec)-> PolynomialVec{
        if self.vec.len()!=v.vec.len(){
            panic!("Mismatched vectors size.")
        }
        let vec=self.vec.iter().zip(v.vec.iter())
        .map(|(b, v)| Polynomial::new(self.n,self.q,&b.coeffs + &v.coeffs))
        .collect::<Vec<Polynomial>>();
        PolynomialVec{n:self.n,q:self.q,vec:vec}

    }
}
/// Rust build-in trait to perform the multiplication of two polynomial_vec 
///instances (scalar product of vector of polynomials)
///
///
////// Arguments:
///
/// V:&PolynomialVec: Reference to a PolynomialVec instance.
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p1_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p2_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p3=&p1*&p2;

impl Mul<&PolynomialVec> for &PolynomialVec{
    type Output=Polynomial;
    fn mul(self, v:&PolynomialVec)-> Polynomial{
        if self.vec.len()!=v.vec.len(){
            panic!("Mismatched vectors size.")
        }
        let vec=self.vec.iter().zip(v.vec.iter()).map(|(x, y)| x * y).collect::<Vec<Polynomial>>();
        let mut acc_pol=Polynomial::new(self.n,self.q,Vec1d(vec![0]));
        for i in 0..vec.len(){
            acc_pol=&acc_pol+&vec[i];
        }
        acc_pol

    }
}
/// Rust build-in trait to perform the multiplication of a PolynomialVec 
///and Vec1d instances (scalar product of vectors in Zq^n and Rq^n)
///
///
////// Arguments:
///
/// V:&PolynomialVec: Reference to a PolynomialVec instance.
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p1_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p2_v=Vec1d(vec![1,2,3]);
/// let p3=&p1*&p2;

impl Mul<&Vec1d> for &PolynomialVec{
    type Output=Polynomial;
    fn mul(self, v:&Vec1d)-> Polynomial{
        if self.vec.len()!=v.0.len(){
            panic!("Mismatched vectors size.")
        }
        let vec=self.vec.iter().zip(v.0.iter()).map(|(x, y)| x * *y).collect::<Vec<Polynomial>>();
        let mut acc_pol=Polynomial::new(self.n,self.q,Vec1d(vec![0]));
        for i in 0..vec.len(){
            acc_pol=&acc_pol+&vec[i];
        }
        acc_pol

    }
}
/// Rust build-in trait to perform the multiplication of a Vec1d
///and PolynomialVec instances (scalar product of vectors in Zq^n and Rq^n)
///
///
////// Arguments:
///
/// V:&PolynomialVec: Reference to a PolynomialVec instance.
///
///
/// Return:
///
/// Polynomial instance.
///
/// Basic usage:
///
/// let p1_v=PolynomialVec::new(128,4099,pol_vec1);
/// let p2_v=Vec1d(vec![1,2,3]);
/// let p3=&p2*&p1;

impl Mul<&PolynomialVec> for &Vec1d{
    type Output=Polynomial;
    fn mul(self, v:&PolynomialVec)-> Polynomial{
        if v.vec.len()!=self.0.len(){
            panic!("Mismatched vectors size.")
        }
        let vec=v.vec.iter().zip(self.0.iter()).map(|(x, y)| x * *y).collect::<Vec<Polynomial>>();
        let mut acc_pol=Polynomial::new(v.n,v.q,Vec1d(vec![0]));
        for i in 0..vec.len(){
            acc_pol=&acc_pol+&vec[i];
        }
        acc_pol
    }
}
/// Function to get a uniformely distributed 
/// Polynomial in Zq[x]/<X^n+1>  
///
/// Arguments:
///
/// * `n`: Length of the field
/// * `q`: Q that generates the field
///
///
///
/// Return:
///
///  Polynomial with the polynomial coefficients sampled

pub fn normal_pol(n:i64,q:i64)->Polynomial{
    Polynomial::new(n,q,normal_vec1d(n as usize,q))
}


/// Function to get a guassian Polynomial in Zq[x]/<X^n+1>  accoding to a gaussian distribution by means of the knuth-yao algorithm.
///
/// Arguments:
///
/// * `n`: Length of the field
/// * `q`: Q that generates the field
///* `mu`: Mean of the gaussian distribution.
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// * `precision`: Number of decimals for the gaussian probabilities.

///
///
/// Return:
///
///  Polynomial with the polynomial coefficients sampled
pub fn gauss_pol(n:i64, q:i64,mu:f64,sigma:f64,t:i64, precision:usize)->Polynomial{
    Polynomial::new(n,q,gauss_vec1d(mu,sigma,t,precision,n as usize))
}
/// Function to get the binary representation of an integer
///
/// Arguments:
///
/// * `value`: Integer from which we want to get the binary representation,
///
///
///
/// Return:
///
///  Vec1d instance whose elements are the elements of the binary representation.
fn int_to_binary(value:i64)->Vec1d{
    let value_str=format!("{:b}", value);
    Vec1d(value_str.chars().map(|e| e.to_digit(10).unwrap() as i64).collect::<Vec<i64>>())
}
/// Definition of BitDecomp trait. Given a polynomial, it returns a polynomial vec with each binary polynomial
/// representation:https://crypto.stackexchange.com/questions/39384/bit-decomposing-a-polynomial-in-bgv-cryptosystem


pub trait BitDecomp{
    fn bitdecomp(&self)-> PolynomialVec ;
}
/// Implementation of BitDecomp for polynomial. It computes the mod q of each coefficient.
///
/// Basic usage:
///
/// let p=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let p=p.bitdecomp();
///
/// Return:
///
/// PolynomialVec instance.
impl BitDecomp for Polynomial{
    fn bitdecomp(&self)->PolynomialVec{
        let l=(self.q as f32).log2().ceil() as usize;
        let coeff_mat=Vec2d{
            N:self.n as usize,
            M:l,
            matrix :self.coeffs.0.clone().into_iter().map(|e| int_to_binary(e).pad(0,l,1).0).collect::<Vec<Vec<i64>>>()
        }.transpose();
        let pol_vec=coeff_mat.matrix.into_iter().rev().map(|e| Polynomial::new(self.n,self.q,Vec1d(e))).collect::<Vec<Polynomial>>();
        PolynomialVec::new(self.n,self.q,pol_vec)

}
}


















