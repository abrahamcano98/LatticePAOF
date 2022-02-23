/* This file contains the necessary traits to performs polynomial operations over the ring Zq[x]/<x^n+1>*/
#![allow(non_snake_case)]

use std::ops::{Mul,Add,Sub};
use crate::linearalgebra::{Vec1d,Modq as Modqla,Pad,circulant_mat,MattoVec,normal_vec1d,gauss_vec1d};

///Polynomial is the struct that represents the polynomials over the finite field Zq[x]/<x^n+1>

#[derive(Debug,Clone)]
pub struct Polynomial{
    ///n: Length of the field
    ///q: q that generates the field
    ///coeffs: coeffs of a polynomial.

    pub n:i64,
    pub q:i64,
    pub coeffs:Vec1d,}

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
        Polynomial::new(self.n,self.q,(&self.coeffs+ &(v.coeffs)).mod_q(self.q))
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
        Polynomial::new(self.n,self.q,(&self.coeffs- &v.coeffs).mod_q(self.q))
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
        Polynomial::new(v.n,v.q,(self*&v.coeffs).mod_q(v.q))
        
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
        Polynomial::new(self.n,self.q,(v*&self.coeffs).mod_q(self.q))
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
/// let p1=Polynomial::new(128,4099,Vec1d(vec![1,0,0]));
/// let p2=Polynomial::new(128,4099,Vec1d(vec![1,1,0,3,4,6]));
/// let p3=&p1*&p2

impl Mul<&Polynomial> for &Polynomial{
    type Output=Polynomial;
    fn mul(self, v:&Polynomial)->Polynomial{
        Polynomial::new(self.n,self.q,(&v.coeffs*&circulant_mat(&self.coeffs)).mattovec().mod_q(self.q))
    }
}
/// Function to sample a uniformely distributed polynomial over Zq[x]/<X^n+1>
/// Arguments:
///
/// * `n`: Length of the field
/// * `q`: Q that generates the field
///
///
/// Return:
///
///  Polynomial with the normal polynomial sampled
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
///  Polynomial with the polynomial sampled
pub fn gauss_pol(n:i64, q:i64,mu:f64,sigma:f64,t:i64, precision:usize)->Polynomial{
    Polynomial::new(n,q,gauss_vec1d(mu,sigma,t,precision,n as usize))
}











