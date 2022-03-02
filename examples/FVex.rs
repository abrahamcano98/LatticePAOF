use PAOF::polynomial::{Polynomial};
use PAOF::linearalgebra::{Vec1d,Modq};
use PAOF::FV::{FV,KeyGen,Encryption,Decryption,Relin};


fn main(){
    let fv=FV{n:16,
        q:12289,
        mu:0.0,
        sigma:3.0,
        t:27,
        precision:120,
        basis:2,};
     
    let key=fv.keygen();

    

    let m0=Polynomial::new(16,12289,Vec1d(vec![0,1,1,1,0,0,0]));
    let m1=Polynomial::new(16,12289,Vec1d(vec![0,1,1,1,0,0,1]));
    println!("Plaintext m0:{:?}",m0.coeffs.0);
    println!("Plaintext m1:{:?}",m1.coeffs.0);

    let m2=Polynomial::new(16,12289,(&m0+&m1).coeffs.mod_q(2));
    let m3=Polynomial::new(16,12289,(&m0*&m1).coeffs.mod_q(2));
    println!("Plaintext m2: m0+m1{:?}",m2.coeffs.0);
    println!("Plaintext m3: m0*m1{:?}",m3.coeffs.0);
    
    let c0=fv.encryption(&key.pk,&m0);
    let c1=fv.encryption(&key.pk,&m1);
    let c2=&c0+&c1;
    let mdec=fv.decryption(&key.sk,&c2);
    println!("Mdec(c0+c1):{:?}",mdec.coeffs.0);

    

    let c3=&c0*&c1;
    let c3p=c3.relin(&key.rlk);
    let mdec=fv.decryption(&key.sk,&c3p);

    println!("Mdec(c0*c1):{:?}",mdec.coeffs.0);



}
