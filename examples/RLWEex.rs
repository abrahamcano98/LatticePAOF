use PAOF::linearalgebra::{Vec1d,Modq};
use PAOF::RLWE::{RLWE as RLWE_prot,KeyGen,Encryption,Decryption};
use PAOF::polynomial::{Polynomial};

fn main(){
    println!("Ring learning with errors encryption scheme.");
        let rlwe=RLWE_prot{    n:32,
            q:12289,
            mu:0.0,
            sigma:3.0,
            t:27,
            precision:120,};
         
        let key=rlwe.keygen();
        let m0=Polynomial::new(32,12289,Vec1d(vec![1,0,1,1,1,0,1,1,1]));
        let m1=Polynomial::new(32,12289,Vec1d(vec![0,1,0,1,1,1,0,1,1]));
        let m2=Polynomial::new(32,12289,(&m0+&m1).coeffs.mod_q(2));
        println!("Plaintext m0:{:?}",m0);
        println!("Plaintext m1:{:?}",m1);
        println!("Plaintex m0+m1:{:?}",m2);
        let c0=rlwe.encryption(&key.pk,&m0);
        let c1=rlwe.encryption(&key.pk,&m1);
        let c2=&c0+&c1;
        let mdec=rlwe.decryption(&key.sk,&c2);
        println!("Decryption of c0+c1 {:?}",mdec);
}