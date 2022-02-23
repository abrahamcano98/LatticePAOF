#![allow(non_snake_case)]

use PAOF::ZKPK::{ZKPK,Prover,Verifier,Proverst,Verifierst};
use PAOF::linearalgebra::{normal_vec2d};
fn main(){
    println!("Amortized zero knowdledge proof of knowledge based on SIS hardness.");
    let zkpk=ZKPK{
        lambda:128,
        q:4099,
        r:5,
        l:4,
        v:4,
        sigma:3.0,
        n:2,
    };
    println!("Preimages we wish to prove we have knowdlege about");
    let S=normal_vec2d(zkpk.v,zkpk.l,zkpk.q);
    println!("S: {:?}",S);
    let prover=Proverst::new(&zkpk,&S);
    let verifier=Verifierst::new(&prover);
    let (W,Y)=prover.compute_w(&zkpk);
    let C=verifier.compute_c(&zkpk);
    let (Z,abort)=prover.compute_z(&C,&Y,&zkpk,2.0);
    if abort==true{
        panic!("Proof aborted. Please run again");
    }
    let result=verifier.verify(&Z,&C,&W,&zkpk); 
    match result{
        true=>println!("Succesful proof. Prover has the knowdledge of the secret S"),
        false=>println!("Succesful proof. Prover does not have the knowdledge of the secret S"),
    }
        
}