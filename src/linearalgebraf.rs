use std::ops::{Mul,Add,Sub};
use crate::linearalgebra::{Vec1d,Vec2d,normal_vec1d};

#[derive(Debug,Clone)]
pub struct Vec1df(pub Vec<f64>);

#[derive(Debug,Clone)]
pub struct Vec2df{
    pub N: usize,
    pub M: usize,
    pub matrix: Vec<Vec<f64>>,
}



impl Mul<&Vec1df> for f64 {
    type Output = Vec1df;
    fn mul(self, v: &Vec1df) -> Vec1df {
        Vec1df(v.0.iter().map(|v| v * self).collect())
    }
}
impl Mul<f64> for &Vec1df {
    type Output = Vec1df;
    fn mul(self, v: f64) -> Vec1df {
        Vec1df(self.0.iter().map(|a| a * v).collect())
    }
}

impl Add<&Vec1df> for &Vec1df{
    type Output=Vec1df;
    fn add(self, v:&Vec1df)-> Vec1df{
        if self.0.len()!=v.0.len(){
            panic!("Mismatched vectors size.")
        }
        Vec1df(self.0.iter().zip(v.0.iter()).map(|(&b, &v)| b + v).collect())
    }
}
impl Sub<&Vec1df> for &Vec1df{
    type Output=Vec1df;
    fn sub(self, v:&Vec1df)-> Vec1df{
        if self.0.len()!=v.0.len(){
            panic!("Mismatched vectors size.")
        }
        Vec1df(self.0.iter().zip(v.0.iter()).map(|(&b, &v)| b - v).collect())
    }
}
impl Mul<&Vec1df> for &Vec1df{
    type Output=f64;
    fn mul(self, v:&Vec1df)-> f64{
        if self.0.len()!=v.0.len(){
            panic!("Mismatched vectors size.")
        }
        self.0.iter().zip(v.0.iter()).map(|(x, y)| x * y).sum::<f64>()
    }
}

impl Mul<&Vec2df> for f64{
    type Output=Vec2df;
    
    fn mul(self, v:&Vec2df)->Vec2df{
    Vec2df{M:v.M,
        N:v.N,
        matrix:(0..v.N)
        .map(|i| {
            (0..v.M)
                .map(|j| v.matrix[i][j]*self).collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>()
       }
    }
}
impl Mul<f64> for &Vec2df{
    type Output=Vec2df;
    
    fn mul(self, v:f64)->Vec2df{
    Vec2df{M:self.M,
        N:self.N,
        matrix:(0..self.N)
        .map(|i| {
            (0..self.M)
                .map(|j| self.matrix[i][j]*v).collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>()
       }
    }

}
impl Add<&Vec2df> for &Vec2df{
    type Output=Vec2df;
    
    fn add(self, v:&Vec2df)->Vec2df{
    if self.M!=v.M || self.N!=v.N{
            panic!("Mismatched matrix size")
    }
    Vec2df{M:self.M,
        N:self.N,
        matrix:(0..self.N)
        .map(|i| {
            (0..self.M)
                .map(|j| self.matrix[i][j]+v.matrix[i][j]).collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>()
       }
    }
}
impl Sub<&Vec2df> for &Vec2df{
    type Output=Vec2df;
    
    fn sub(self, v:&Vec2df)->Vec2df{
    if self.M!=v.M || self.N!=v.N{
            panic!("Mismatched matrix size")
    }
    Vec2df{M:self.M,
        N:self.N,
        matrix:(0..self.N)
        .map(|i| {
            (0..self.M)
                .map(|j| self.matrix[i][j]-v.matrix[i][j]).collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>()
       }
    }
}
impl Mul<&Vec2df> for &Vec2df{
    type Output=Vec2df;
    
    fn mul(self, v:&Vec2df)->Vec2df{
    if self.M!=v.N{
            panic!("Mismatched matrix size")
    }
        Vec2df{
            N:self.N, 
            M:v.M ,
            matrix:
            (0..self.N)
            .map(|i| {
                (0..v.M)
                    .map(|j| (0..self.M).map(|k| self.matrix[i][k] * v.matrix[k][j]).sum())
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<Vec<f64>>>()
        }

}
}
impl Mul<&Vec2df> for &Vec1df{
    type Output=Vec2df;
    
    fn mul(self, v:&Vec2df)->Vec2df{
    &self.vectomat()*v
   
}
}

impl Mul<&Vec1df> for &Vec2df{
    type Output=Vec2df;

    fn mul(self, v:&Vec1df)->Vec2df{
        self*&v.transpose()
    }
}

pub trait Sum<T>{
    fn sum(&self)->T;
}

impl Sum<f64> for Vec1df{
    fn sum(&self)->f64{
        self.0.iter().sum()
    }
}
impl Sum<f64> for Vec2df{
    fn sum(&self)->f64{
        self.matrix.iter().map(|x| -> f64 { x.iter().sum() }).sum()
    }
}
pub trait Norm<U>{
    fn norm(&self)->U;
}

impl Norm<f64> for Vec1df{
    fn norm(&self)->f64{
        (0..self.0.len()).map(|i| (self.0[i]*self.0[i]) as f64).sum::<f64>().sqrt()
    }
}

pub trait Transpose<T>{
    fn transpose(&self)->T;
}

impl Transpose<Vec2df> for Vec1df{
    fn transpose(&self)->Vec2df{
        Vec2df{
            M:1,
            N:self.0.len(),
            matrix:(0..self.0.len())
        .map(|i| {
            (0..1)
                .map(|_| self.0[i]).collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>()
    }

}}
impl Transpose<Vec2df> for Vec2df{
    fn transpose(&self)->Vec2df{
        Vec2df{
            M:self.N,
            N:self.M,
            matrix:(0..self.M)
        .map(|i| {
            (0..self.N)
                .map(|j| self.matrix[j][i]).collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>()
    }

}}

pub trait VectoMat<T>{
    fn vectomat(&self)->T;
}

impl VectoMat<Vec2df> for Vec1df{
    fn vectomat(&self)->Vec2df{
        let mut new_vec:Vec<Vec<f64>>=Vec::new();
        new_vec.push(self.0.clone());
        Vec2df{
            N:1,
            M:self.0.len(),
            matrix:new_vec
        }
    }
}

pub trait MattoVec<T>{
    fn mattovec(&self)->T;
}

impl MattoVec<Vec1df> for Vec2df{
    fn mattovec(&self)->Vec1df{
        if self.N!=1{
            panic!("Cannot be converted any matrix with more than 1 rows to a vector")
        }
        Vec1df(self.matrix[0].clone())
    }
}

pub fn power_iteration(mut A:Vec2d,max_iter:i64)->f64{
    let mut b_k=Vec1d(vec![5;A.M]);
    let mut b_k=Vec1df(b_k.0.iter().map(|i| *i as f64).collect::<Vec<f64>>()).transpose();
    let B=Vec2df{N:A.N,M:A.M,matrix:(0..A.N).map(|i| {(0..A.M).map(|j| A.matrix[i][j] as f64 ).collect::<Vec<f64>>()}).collect::<Vec<Vec<f64>>>()};
    for i in 0..max_iter{
        let b_k1=(&B*(&b_k));
        let b_k1_norm=b_k1.transpose().mattovec().norm();
        b_k = (1.0/b_k1_norm)*&b_k1;
    }
    let ray_num=(&(&b_k.transpose()*&B)*&b_k);
    let ray_div=(&(b_k.transpose())*&b_k);
    println!("result:{}",(ray_num.matrix[0][0]/ray_div.matrix[0][0]).sqrt().abs());
    (ray_num.matrix[0][0]/ray_div.matrix[0][0]).sqrt().abs()
}

