#![allow(non_snake_case)]

use PAOF::ZKPK::{ZKPK,Prover,Verifier,Proverst,Verifierst};
use PAOF::linearalgebra::{normal_vec2d};
fn main(){
    println!("Amortized zero knowdledge proof of knowledge based on SIS hardness.");
    
    let zkpk=ZKPK::new(128,4099,5,4,7,3.0);
    
    println!("Preimages we wish to prove we have knowdlege about");
    
    let S=normal_vec2d(zkpk.v,zkpk.l,2);
    
    println!("S: {:?}",S);
    
    let prover=Proverst::new(&zkpk,&S);
    let verifier=Verifierst::new(prover.A.clone(),prover.T.clone(),prover.B.clone());
    println!("Computing proof...");

    let (W,Y)=prover.compute_w(&zkpk);
    let C=verifier.compute_c(&zkpk);
    let (Z,abort)=prover.compute_z(&C,&Y,zkpk.rho);
    if abort==true{
        panic!("Proof aborted. Please run again");
    }
    let result=verifier.verify(&Z,&C,&W,&zkpk); 
    match result{
        true=>println!("Succesful proof. Prover has the knowdledge of the secret S"),
        false=>println!("Succesful proof. Prover does not have the knowdledge of the secret S"),
    }
        
}
