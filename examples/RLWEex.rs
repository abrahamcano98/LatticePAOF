use PAOF::linearalgebra::{Vec1d,Vec2d,Transpose,Pad,normal_vec2d,gauss_vec2d};
use PAOF::RLWE::{Key,RLWE as RLWE_prot,KeyGen,Encryption,Decryption};
use PAOF::polynomial::Polynomial;

fn main(){
        let rlwe=RLWE_prot{    n:512,
            q:12289,
            mu:0.0,
            sigma:3.0,
            t:27,
            precision:120,};
        let key=rlwe.keygen();
        let m0=Polynomial::new(512,12289,Vec1d(vec![1,0,1,0,0,0,1]));
        let c0=rlwe.encryption(key.pk,m0.clone());
        let mdec=rlwe.decryption(key.sk,c0);
        println!("m0:{:?}",mdec);
        println!("mdec:{:?}",m0);
        }