use crate::linearalgebra::{Vec2d,normal_vec2d,Transpose,Pad,Modq,Sum};
use crate::polynomial::Polynomial;


///CommitKey is the structure that represents the components of the commitment key based on SIS problem.
pub struct CommitKey{
   pub A1:Vec2d,
   pub A2:Vec2d,
}
///Sis commitment is the structure that represents the parameters used by the SIS commitment protocol
pub struct SisCommitment{
    /// *`q`: q that generates the field.
    /// *`Q`: New generator.
    /// *`r`: Rows of commitment keys.
    ///*`n`: Number of elements to be sent.
    q:i64,
    Q:i64,
    r:usize,
    n:usize,
}
/// Definition of sis commitment constructor.
impl SisCommitment{
    pub fn new(q:i64, r:usize,n:usize) -> SisCommitment {
        SisCommitment {
            q: q,
            Q: next_prime(q),
            r:r,
            n:n
        }
    }
}
/// Definition of the trait KeyGen. It generates the commitment key of the protocol.

pub trait KeyGen<T> {
    fn keygen(&self)->T;
}
/// Implementation of trait KeyGen for RLWE. It samples A1 from a normal distribution in Zq^{r x 2*logq Q}
/// It samples A2 from a normal distribution in Zq^{r x n}
impl KeyGen<CommitKey> for SisCommitment{
    fn keygen(&self)->CommitKey{
        let log=(2.0_f64*((self.Q as f64).log(self.q as f64))).ceil() as usize;
        CommitKey{
            A1:normal_vec2d(self.r,log,self.Q),
            A2:normal_vec2d(self.r, self.n,self.Q),
        }
    }
}
/// Definition of the trait Commitment. It commits messages of length n in Zq^n and outputs commitments c in Z_{Q}^r


pub trait Commitment<T> {
    fn commitment(&self,ck:&CommitKey,m:&T)-> (T,T);
}
/// Implementation of commitment for Sis Commitment
///c=(A1*r+A2*m) mod Q
///
/// Arguments:
///
///* `ck`: Commitment key.
/// * `m`: Vec2d with represents the matrix message to be commited.
///
///
/// Return:
///
///  (Commitment, randomness)

impl Commitment<Vec2d> for SisCommitment{
    fn commitment(&self,ck:&CommitKey,mut m: &Vec2d)-> (Vec2d,Vec2d){
        
        let log=(2.0_f64*((self.Q as f64).log(self.q as f64))).ceil() as usize;
        let r=normal_vec2d(m.clone().N, log, self.q).transpose();
        let mp=&m.clone().pad(0,self.n,0);
        let mt=&mp.transpose().mod_q(self.q);
        let A1=&ck.A1;
        let A2=&ck.A2;

        

        (&(A1*&r)+&(A2*mt).mod_q(self.Q),r)
    }
}
/// Definition of the opening of the commitment.

pub trait OpenComm<T>{
    fn open(&self,ck:&CommitKey,r:&T,c:&T,m:&T)->bool;
}
/// Implementation of commitment for Sis Commitment
///c-(A1*r+A2*m)=[0] mod Q
///
/// Arguments:
///
///* `ck`: Commitment key.
/// * `r`: Randomness used in the commitment.
/// *`c`: Commitment computed.
/// *`m`: Matrix with the plain message.
///
/// Return:
///
///  Boolean. True if the prover sends the original message, false otherwise.
impl OpenComm<Vec2d> for SisCommitment{
    fn open(&self,ck:&CommitKey,r:&Vec2d,c:&Vec2d,m:&Vec2d)->bool{
        let A1=&ck.A1;
        let A2=&ck.A2;
        let mp=&m.clone().pad(0,self.n,0);
        let mt=&mp.transpose().mod_q(self.q);
        (&(c-&(A1*r))-&(A2*mt)).mod_q(self.Q).sum()==0
    }
}



/// Function to know if a number is prime according to the Erastothenes method.
///
/// Arguments:
///
/// * `n`: Number to check
///
/// Return:
///
///  Bool
fn is_prime(n:i64)->bool{
    if n<=3{
        return true;
    }
    if n%2==0 || n%3==0{
       return false;
    }
    let mut i:i64=5;
    while i.pow(2) <=n{
        if n%i==0 || n%(i+2)==0{
            return false;
        }
        i=i+6;
    }
    return true;
}
/// Function to obtain the next prime number given one.
///
/// Arguments:
///
/// * `n`: Prime number
///
/// Return:
///
///  i64 which is the next prime.

fn next_prime(mut n:i64)->i64{
loop{
    n=n+1;
    if is_prime(n)==false{
        return n;
    }
}
}