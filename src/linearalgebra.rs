/*This file contains  the minimal functions and traits to perform linear algebra operations over finite fields"*/

use std::ops::{Mul,Add,Sub};
use crate::samplers;


///Vec1d is the struct that represents the vectors over the finite field Z


#[derive(Debug,Clone)]
pub struct Vec1d(pub Vec<i64>);
///Vec2d is the struct that represents a matrix over the finite field Z


#[derive(Debug,Clone)]
pub struct Vec2d{
    ///M: Number of columns
    ///N: Number of rows
    ///Matrix: Matrix

    pub N: usize,
    pub M: usize,
    pub matrix: Vec<Vec<i64>>,
}

/// Rust build-in trait to perform the multiplication of a vector by a scalar. 
///
///
////// Arguments:
///
/// V:&Vec1d: Reference to a Vec1d struct instance

///
/// Return:
///
/// Vec1d instance
///
/// Basic usage:
///let mut v=Vec1d(vec![1,2,3]);
///v5=5_i64*&v;
impl Mul<&Vec1d> for i64 {
    type Output = Vec1d;
    fn mul(self, v: &Vec1d) -> Vec1d {
        Vec1d(v.0.iter().map(|v| v * self).collect())
    }
}
/// Rust build-in trait to perform the multiplication of a vector by a scalar. 
///
///
////// Arguments:
///
/// V:i64: i64 number

///
/// Return:
///
/// Vec1d instance
///
/// Basic usage:
///let v=Vec1d(vec![1,2,3]);
///let v5=&v*5_i65; 
impl Mul<i64> for &Vec1d {
    type Output = Vec1d;
    fn mul(self, v: i64) -> Vec1d {
        Vec1d(self.0.iter().map(|a| a * v).collect())
    }
}
/// Rust build-in trait to perform the addition of two Vec1d instances 
///
///
////// Arguments:
///
/// V:&Vec1d: Reference to a Vec1d instance.

///
/// Return:
///
/// Vec1d instance.
///
/// Basic usage:
///
///let v1=Vec1d(vec![1,2,3]);
///let v2=Vec1d(vec![4,5,6]):
///let v3=&v1+&v2

impl Add<&Vec1d> for &Vec1d{
    type Output=Vec1d;
    fn add(self, v:&Vec1d)-> Vec1d{
        if self.0.len()!=v.0.len(){
            panic!("Mismatched vectors size.")
        }
        Vec1d(self.0.iter().zip(v.0.iter()).map(|(&b, &v)| b + v).collect())
    }
}

/// Rust build-in trait to perform the substraction of two Vec1d instances 
///
///
////// Arguments:
///
/// V:&Vec1d: Reference to a Vec1d instance.

///
/// Return:
///
/// Vec1d instance.
///
/// Basic usage:
///
///let v1=Vec1d(vec![1,2,3]);
///let v2=Vec1d(vec![4,5,6]):
///let v3=&v1-&v2


impl Sub<&Vec1d> for &Vec1d{
    type Output=Vec1d;
    fn sub(self, v:&Vec1d)-> Vec1d{
        if self.0.len()!=v.0.len(){
            panic!("Mismatched vectors size.")
        }
        Vec1d(self.0.iter().zip(v.0.iter()).map(|(&b, &v)| b - v).collect())
    }
}
/// Rust build-in trait to perform the multiplication of two Vec1d instances (Scalar product: <V,W>=a) 
///
///
////// Arguments:
///
/// V:&Vec1d: Reference to a Vec1d instance.

///
/// Return:
///
/// Vec1d instance.
///
/// Basic usage:
///
///let v1=Vec1d(vec![1,2,3]);
///let v2=Vec1d(vec![4,5,6]):
///let v3:i64=&v1*&v2
impl Mul<&Vec1d> for &Vec1d{
    type Output=i64;
    fn mul(self, v:&Vec1d)-> i64{
        if self.0.len()!=v.0.len(){
            panic!("Mismatched vectors size.")
        }
        self.0.iter().zip(v.0.iter()).map(|(x, y)| x * y).sum::<i64>()
    }
}
/// Rust build-in trait to perform the multiplication of a matrix by a scalar. 
///
///
////// Arguments:
///
/// V:&Vec2d Reference to a Vec2d instance

///
/// Return:
///
/// Vec2 instance
///
/// Basic usage:
///let v=Vec2d(vec![vec![1,2],vec![3,4]]);
///let v5=5_i64*v; 

impl Mul<&Vec2d> for i64{
    type Output=Vec2d;
    
    fn mul(self, v:&Vec2d)->Vec2d{
    Vec2d{M:v.M,
        N:v.N,
        matrix:(0..v.N)
        .map(|i| {
            (0..v.M)
                .map(|j| v.matrix[i][j]*self).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
       }
    }
}

/// Rust build-in trait to perform the multiplication of a matrix by a scalar. 
///
///
////// Arguments:
///
/// V:i64 i64 number

///
/// Return:
///
/// Vec2 instance
///
/// Basic usage:
///let v=Vec2d(vec![vec![1,2],vec![3,4]]);
///let v5=v*5_i64; 

impl Mul<i64> for &Vec2d{
    type Output=Vec2d;
    
    fn mul(self, v:i64)->Vec2d{
    Vec2d{M:self.M,
        N:self.N,
        matrix:(0..self.N)
        .map(|i| {
            (0..self.M)
                .map(|j| self.matrix[i][j]*v).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
       }
    }

}
/// Rust build-in trait to perform the addition of two Vec2d instances 
///
///
////// Arguments:
///
/// V:&Vec2d: Reference to a Vec2d instance.

///
/// Return:
///
/// Vec2d instance.
///
/// Basic usage:
///
///let v1=Vec2d(vec![vec![1,2],vec![3,4]]);
///let v2=Vec2d(vec![vec![5,6],vec![7,8]]);
///let v3=&v1+&v2



impl Add<&Vec2d> for &Vec2d{
    type Output=Vec2d;
    
    fn add(self, v:&Vec2d)->Vec2d{
    if self.M!=v.M || self.N!=v.N{
            panic!("Mismatched matrix size")
    }
    Vec2d{M:self.M,
        N:self.N,
        matrix:(0..self.N)
        .map(|i| {
            (0..self.M)
                .map(|j| self.matrix[i][j]+v.matrix[i][j]).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
       }
    }
}
/// Rust build-in trait to perform the substraction of two Vec2d instances 
///
///
////// Arguments:
///
/// V:&Vec2d: Reference to a Vec2d instance.

///
/// Return:
///
/// Vec2d instance.
///
/// Basic usage:
///
///let v1=Vec2d(vec![vec![1,2],vec![3,4]]);
///let v2=Vec2d(vec![vec![5,6],vec![7,8]]);
///let v3=&v1-&v2
impl Sub<&Vec2d> for &Vec2d{
    type Output=Vec2d;
    
    fn sub(self, v:&Vec2d)->Vec2d{
    if self.M!=v.M || self.N!=v.N{
            panic!("Mismatched matrix size")
    }
    Vec2d{M:self.M,
        N:self.N,
        matrix:(0..self.N)
        .map(|i| {
            (0..self.M)
                .map(|j| self.matrix[i][j]-v.matrix[i][j]).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
       }
    }
}
/// Rust build-in trait to perform the substraction of two Vec2d instances 
///
///
////// Arguments:
///
/// V:&Vec2d: Reference to a Vec2d instance.

///
/// Return:
///
/// Vec2d instance.
///
/// Basic usage:
///
///let v1=Vec2d(vec![vec![1,2],vec![3,4]]);
///let v2=Vec2d(vec![vec![5,6],vec![7,8]]);
///let v3=&v1*&v2
impl Mul<&Vec2d> for &Vec2d{
    type Output=Vec2d;
    
    fn mul(self, v:&Vec2d)->Vec2d{
    if self.M!=v.N{
            panic!("Mismatched matrix size")
    }
        Vec2d{
            N:self.N, 
            M:v.M ,
            matrix:
            (0..self.N)
            .map(|i| {
                (0..v.M)
                    .map(|j| (0..self.M).map(|k| self.matrix[i][k] * v.matrix[k][j]).sum())
                    .collect::<Vec<i64>>()
            })
            .collect::<Vec<Vec<i64>>>()
        }

}
}
/// Rust build-in trait to perform the multiplication of a Vec2d instance by a Vec1d vector 
///
///
////// Arguments:
///
/// V:&Vec2d: Reference to a Vec2d instance.

///
/// Return:
///
/// Vec2d instance.
///
/// Basic usage:
///
///let v1=Vec1d(vec![1,2,3]);
///let v2=Vec2d(vec![vec![1,2],vec![3,4],vec![5,6]]);
///let v3=&v1*&v2
impl Mul<&Vec2d> for &Vec1d{
    type Output=Vec2d;
    
    fn mul(self, v:&Vec2d)->Vec2d{
    &self.vectomat()*v
   
}
}
/// Rust build-in trait to perform the multiplication of a Vec2d instance by a Vec1d vector 
///
///
////// Arguments:
///
/// V:&Vec1d: Reference to a Vec1d instance.

///
/// Return:
///
/// Vec2d instance.
///
/// Basic usage:
///
///let v1=Vec2d(vec![vec![1,2],vec![3,4],vec![5,6]]);
///let v2=Vec1d(vec![1,2])
///let v3=&v1*&v2;
impl Mul<&Vec1d> for &Vec2d{
    type Output=Vec2d;

    fn mul(self, v:&Vec1d)->Vec2d{
        self*&v.transpose()
    }
}
///Trait sum to return the sum of all elements given a Vec1d instance or Vec2d instance.

pub trait Sum<T>{
    fn sum(&self)->T;
}
///Implementation of trait sum for Vec1d
///
/// Basic usage:
///
///let v1=Vec1d(vec![1,2,3]);
///let sum:i64=v1.sum();
impl Sum<i64> for Vec1d{
    fn sum(&self)->i64{
        self.0.iter().sum()
    }
}
///Implementation of trait sum for Vec2d
///
/// Basic usage:
///
///let v1=Vec2d(vec![vec![1,2],vec![3,4],vec![5,6]]);
///let sum:i64=v1.sum();
impl Sum<i64> for Vec2d{
    fn sum(&self)->i64{
        self.matrix.iter().map(|x| -> i64 { x.iter().sum() }).sum()
    }
}
///Trait Modq to convert samples in the Ring Z to the ring Zq


pub trait Modq<T>{
    fn mod_q(&self,q:T)->Self;
}
///Implementation of trait Modq for usize
///
///Arguments
///
///q:usize Ring modulus q
///
/// Basic usage:
///
///let a=5;
///let q=2;
/// let aq=a.mod_q(q);
impl Modq<usize> for usize{
    fn mod_q(&self, q:usize)->usize{
        
        ((self%q)+q)%q 
    }
}
///Implementation of trait Modq for Vec1d
///
///Arguments
///
///q:usize Ring modulus q
///
/// Basic usage:
///
///let a=Vec1d(vec![1,2,3]);
///let q=2;
/// let aq=a.mod_q(q);
impl Modq<i64> for Vec1d{
    fn mod_q(&self,q:i64)->Vec1d{
        Vec1d((0..self.0.len()).map(|i| ((self.0[i] % q) + q) % q).collect::<Vec<i64>>())
    }
}
///Implementation of trait Modq for Vec2d
///
///Arguments
///
///q:i64 Ring modulus q
///
/// Basic usage:
///
///let a=Vec1d(vec![vec![1,2,3],vec![4,5,6]]);
///let q=2;
/// let aq=a.mod_q(q);
impl Modq<i64> for Vec2d{
    fn mod_q(&self,q:i64)->Vec2d{
        Vec2d{N:self.N,
            M:self.M,
            matrix: (0..self.N)
            .map(|i| {
                (0..self.M)
                    .map(|j| ((self.matrix[i][j] % q) + q) % q).collect::<Vec<i64>>()
            })
            .collect::<Vec<Vec<i64>>>()}
    }
}

/// Definition of trait push, which is used to add elements to the structures.
pub trait Push<T>{
    fn push(self, v:T)->Self;
}
///Implementation of trait Push for Vec1d
///
///Arguments
///
///v:i64 Value to be added
///
/// Basic usage:
///
///let a=Vec1d(vec![1,2,3]);
///let v=2;
/// let aq=a.push(v);
impl Push<i64> for Vec1d{
    fn push(mut self,v:i64)->Vec1d{
        self.0.push(v);
        self
        
    }
}
/// Definition of trait push, which is used to add elements in a certain position to the structures.

pub trait Insert<T,U>{
    fn insert(self, i:T,v:U)->Self;
}
///Implementation of trait Push for Vec1d
///
///Arguments
///
///i:usize Position to be added
///v:i64 Value to be added
///
/// Basic usage:
///
///let a=Vec1d(vec![1,2,3]);
///let v=2;
///let i=2;
/// let aq=a.insert(i,v);
impl Insert<usize,i64> for Vec1d{
    fn insert(mut self,i:usize,v:i64)->Self{
        (self.0).insert(i,v);
        self 
    }
}
 
impl Insert<usize,i64> for &mut Vec1d{
    fn insert(self,i:usize,v:i64)->Self{
        (self.0).insert(i,v);
        self 
    }
}
/// Definition of trait pad, which is used to pad vec1d and vec2d structures. 

pub trait Pad<T,U>{
    fn pad(self,v:T,c:U,d:T)->Self;
}
///Implementation of trait Pad for Vec1d
///
///Arguments
///
///
///v:i64 Value to be padded.
///c:usize Total length of the desired vector
///d: i64 rientation of the padding. d=0 to pad through to the right side and d=1 to the left
///
/// Basic usage:
///
///let a=Vec1d(vec![1,2,3]);
/// let ap=a.pad(0,5,0);

impl Pad<i64,usize> for Vec1d{
    fn pad(mut self, v:i64,c:usize ,d:i64)->Vec1d{
        let len=self.0.len();
        if d==0{
            for i in len..c{
                self.0.insert(i,v);
            }
        }
        else if d==1{
            for i in 0..(c-len){
                self.0.insert(i,v);
            }
        }
        else{
            panic!("Please, set either d=0 to pad on the right side or d=1 to pad on the left side")
        }
        self
        
    }
}
///Implementation of trait Pad for Vec1d. It pads each row vector in a similar way to the Vec1d implementation.
///
///Arguments
///
///
///v:i64 Value to be padded.
///c:usize Total length of the desired vector
///d: i64 rientation of the padding. d=0 to pad through to the right side and d=1 to the left
///
/// Basic usage:
///
///let a=Vec2d(vec![vec![1,2],vec![3,4]]);
/// let ap=a.pad(0,5,0);

impl Pad<i64,usize> for Vec2d{
    fn pad(self,v:i64, c:usize, d:i64)->Vec2d{
    let mut matrix=Vec::new();
    for i in 0..self.N{
        matrix.push(Vec1d(self.matrix[i].clone()).pad(v,c,d).0)
    }
    Vec2d{N:self.N,M:c,matrix:matrix}
}
}
/// Definition of trait Norm. It is used to compute the euclidean norm of a vector.
pub trait Norm<U>{
    fn norm(&self)->U;
}
///Implementation of trait Norm for Vec1d. It pads each row vector in a similar way to the Vec1d implementation.
///
///return: f64.
///
/// Basic usage:
///
///let a=Vec1d(vec![1,2,3]);
/// let norm=a.norm();
impl Norm<f64> for Vec1d{
    fn norm(&self)->f64{
        (0..self.0.len()).map(|i| (self.0[i]*self.0[i]) as f64).sum::<f64>().sqrt()
    }
}
/// Definition of trait Flatten. It is used to express matrixes as a vectors.
pub trait Flatten<T>{
    fn flatten(&self)->T;
}
///Implementation of trait Flatten for Vec2d. 
///
///
///return: Vec1d.
///
///
/// Basic usage:
///
///let a=Vec2d(vec![vec![1,2],vec![3,4]]);
/// let norm=a.flatten();
impl Flatten<Vec1d> for Vec2d{
    fn flatten(&self)->Vec1d{
        let flatten_array: Vec<i64> = self.matrix
        .iter()
        .flat_map(|array| array.iter())
        .cloned()
        .collect();
        Vec1d(flatten_array)
}
}
/// Definition of trait transpose, which is used to transpose both vectors and matrixes

pub trait Transpose<T>{
    fn transpose(&self)->T;
}
///Implementation of trait Transpose for Vec2d. 
///
///
///return: Vec2d.
///
///
/// Basic usage:
///
///let a=Vec1d(vec![1,2,3]);
/// let at=a.transpose();

impl Transpose<Vec2d> for Vec1d{
    fn transpose(&self)->Vec2d{
        Vec2d{
            M:1,
            N:self.0.len(),
            matrix:(0..self.0.len())
        .map(|i| {
            (0..1)
                .map(|_| self.0[i]).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
    }

}}
///Implementation of trait Transpose for Vec2d. 
///
///
///return: Vec2d.
///
///
/// Basic usage:
///
///let a=Vec2d(vec![vec![1,2],vec![3,4]]);
/// let at=a.transpose();
impl Transpose<Vec2d> for Vec2d{
    fn transpose(&self)->Vec2d{
        Vec2d{
            M:self.N,
            N:self.M,
            matrix:(0..self.M)
        .map(|i| {
            (0..self.N)
                .map(|j| self.matrix[j][i]).collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>()
    }

}}
/// Definition of trait VecToMat. It is used to convert Vec1d structures to Vec2d.
pub trait VectoMat<T>{
    fn vectomat(&self)->T;
}
///Implementation of trait VectoMat for Vec1d. 
///
///
///return: Vec2d.
///
///
/// Basic usage:
///
///let a=Vec1d(vec![1,2,3]);
/// let norm=a.vectomat();
impl VectoMat<Vec2d> for Vec1d{
    fn vectomat(&self)->Vec2d{
        let mut new_vec:Vec<Vec<i64>>=Vec::new();
        new_vec.push(self.0.clone());
        Vec2d{
            N:1,
            M:self.0.len(),
            matrix:new_vec
        }
    }
}
/// Definition of trait MatToVec. It is used to convert Vec2d structures to Vec1d.
pub trait MattoVec<T>{
    fn mattovec(&self)->T;
}
///Implementation of trait MatToVec for Vec12d. 
///
///
///return: Vec1d.
///
///
/// Basic usage:
///
///let a=Vec2d(vec![vec![1,2,3],vec![4,5,6]]);
/// let norm=a.vectomat();
impl MattoVec<Vec1d> for Vec2d{
    fn mattovec(&self)->Vec1d{
        if self.N!=1{
            panic!("Cannot be converted any matrix with more than 1 rows to a vector")
        }
        Vec1d(self.matrix[0].clone())
    }
}
/// Function to concatenate two Vec1d vectors.
///
///
/// Arguments:
///
/// * `a`: Vec1d
/// * 'b': Vec1d

///
/// Return:
///
///  Vec1d which contains both vectors concatenated.
pub fn concatenate_vec1d(mut a: Vec1d, b:Vec1d)->Vec1d{
    (a.0).extend(&(b.0));
    a
}
/// Function to concatenate two compute the circulant matrix created by a vector.
///
///
/// Arguments:
///
/// * `x`: Vec1d

///
/// Return:
///
///  Vec2d which contains the matrix generated by x

pub fn circulant_mat(x:&Vec1d)->Vec2d{
    let len=x.0.len();
    Vec2d{N:len,
    M:len,
    matrix: (0..len)
    .map(|i| {
        (0..len)
            .map(|j| sign(i,j)*x.0[((j as i64-i as i64) as usize).mod_q(len)]).collect::<Vec<i64>>()
    })
    .collect::<Vec<Vec<i64>>>()}
}


fn sign(i:usize, j:usize)->i64{
    if i>j{
        -1
    }
    else{
        1
    }
}

/// Function to sample a uniformely distributed vector over Zq.
///
/// Arguments:
///
/// * `size`: Size of the vector
/// * `max_value`: q which generates the field

///
/// Return:
///
///  Vec1d with the normal vector sampled
pub fn normal_vec1d(size:usize, max_value:i64)->Vec1d{
    Vec1d(samplers::normal_vec(size,max_value))
}
/// Function to sample a uniformely distributed matrix over Zq.
///
/// Arguments:
///
/// * `height`: Rows of the matrix
/// * `width`: Columns of the matrix
/// * `max_value`: q which generates the field

///
/// Return:
///
///  Vec2d with the normal matrix sampled

pub fn normal_vec2d(height:usize, width:usize, max_value:i64)->Vec2d{
    Vec2d{N:height,M:width,matrix:samplers::normal_mat(height,width,max_value)}
}
/// Function to sample a discrete gaussian vector over Zq.
///
/// Arguments:
///
/// * `mu`: Mean of the gaussian.
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// * `precision`: Gaussian probabilities decimals. The greater the number the more approximate the gaussian.
/// * `length`: Length of the vector
///
///
/// Return:
///
///  Vec1d with the discrete gaussian vector sampled
pub fn gauss_vec1d(mu:f64,sigma:f64,t:i64, precision:usize, length:usize)->Vec1d{
    Vec1d(samplers::gauss_vec(mu,sigma,t,precision,length))
}
/// Function to sample a discrete gaussian matrix over Zq.
///
/// Arguments:
///
/// * `mu`: Mean of the gaussian.
/// * `sigma`: Standard deviation of the gaussian.
/// * `t`: Max value to be sampled according the sigma: (-t*sigma..t*sigma)
/// * `precision`: Gaussian probabilities decimals. The greater the number the more approximate the gaussian.
/// * `length`: Rows of the matrix
/// * `width`: Columns of the matrix
///
///
/// Return:
///
///  Vec2 with the discrete gaussian matrix sampled
pub fn gauss_vec2d(mu:f64,sigma:f64,t:i64, precision:usize, length:usize,width:usize)->Vec2d{
    Vec2d{N:length,M:width,matrix:samplers::gauss_mat(mu,sigma,t,precision,length,width)}
}

