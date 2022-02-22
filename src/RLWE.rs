use crate::linearalgebra::Vec1d;
use crate::polynomial::{Polynomial,normal_pol,gauss_pol};

///RLWE is the struct that represents the parameters employed to encrypt a message over Zq[x]/<x^n+1>


pub struct RLWE{
    ///n: Length of the field
    ///q: Q that generates the field
    ///mu: Mean of the error distribution
    ///sigma: Standard deviation of the error distribution
    ///t:Max value to be sampled according the sigma: (-t*sigma..t*sigma)
    ///precision:
    pub n:i64,
    pub q:i64,
    pub mu:f64,
    pub sigma:f64,
    pub t:i64,
    pub precision:usize,
}
/// Key structure to store the keys employed by the cipher
pub struct Key{
    pub sk:Polynomial,
    pub pk:(Polynomial,Polynomial)
}
/// Polynomial structure to store the components employed by the cipher.
pub struct Ciphertext{
    pub u:Polynomial,
    pub v:Polynomial,
}
/// Definition of the trait KeyGen. It generates both sk and pk of the protocol
pub trait KeyGen<T> {
    fn keygen(&self)-> T;
}
/// Implementation of trait KeyGen for RLWE. It samples A from Zq[x]/<x^n+1> and sk,e from the error distribution.
/// The sk is set to be the sk, while the pk is equal to a*sk+e

impl KeyGen<Key> for RLWE{
    fn keygen(&self)->Key{
        let a=normal_pol(self.n,self.q);
        let sk=gauss_pol(self.n, self.q,self.mu,self.sigma,self.t, self.precision);
        let e=gauss_pol(self.n, self.q,self.mu,self.sigma,self.t, self.precision);
        Key{sk:sk.clone(),pk:(a.clone(),&(&a*&sk)+&e)}
    }
}
/// Definition of the trait Encryption. It generates a ciphertext given a plaintext message m and a pk.
pub trait Encryption<T,U> {
    fn encryption(&self,pk:(Polynomial,Polynomial),m:T)-> U;
}
/// Implementation of Encryption for RLWE.
///u=a*r+e1
///v=b*r+e2+[q/2]*m
///
/// Arguments:
///
///* `pk`: A tuple with the components of the pk.
/// * `m`: The plaintext message.
///
///
/// Return:
///
///  Ciphertext structure.
impl Encryption<Polynomial,Ciphertext> for RLWE{
    fn encryption(&self,pk:(Polynomial,Polynomial),m:Polynomial)->Ciphertext{
        let a=pk.0;
        let b=pk.1;
        let r=gauss_pol(self.n, self.q,self.mu,self.sigma,self.t, self.precision);
        let e1=gauss_pol(self.n, self.q,self.mu,self.sigma,self.t, self.precision);
        let e2=gauss_pol(self.n, self.q,self.mu,self.sigma,self.t, self.precision);
        let qhr=(self.q as f64/2.0).ceil() as i64;
       
        let u=&(&a*&r)+&e1; 
        let v=&(&(&b*&r)+&e2)+&(qhr*&m);
        
       
        Ciphertext{u:u,v:v}
    }
}
/// Definition of the trait Decryption. It generates a plaintext given a ciphertext and a secret key.
pub trait Decryption<T,U> {
    fn decryption(&self,sk:Polynomial,c:T)-> U;
}
/// Implementation of Decryption for RLWE.
///
/// mdec=v-u*sk
/// for each mdec round to either 0 or [q/2] whichever is closer.
/// Arguments:
///
///* `sk`: Polynomial with the sk generated by keygen.
/// * `c`: Ciphertext structure
///
///
/// Return:
///
///  Plaintext message.

impl Decryption<Ciphertext,Polynomial> for RLWE{
    fn decryption(&self,sk:Polynomial,c:Ciphertext)->Polynomial{
        let u=c.u;
        let v=c.v;
        let m=&v-&(&u*&sk);
        let qhr=(self.q as f64/2.0).ceil();
        let m_cen=(0..self.n).map(|i| {
            if m.coeffs.0[i as usize]<qhr as i64{
                m.coeffs.0[i as usize]
        }else{
            m.coeffs.0[i as usize]-self.q
        }}).collect::<Vec<i64>>();
        Polynomial::new(self.n,self.q,Vec1d((0..self.n).map(|i| (m_cen[i as usize] as f64/(self.q as f64/4.0)).abs() as i64).collect::<Vec<i64>>()))
    }
}
