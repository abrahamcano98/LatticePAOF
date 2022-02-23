use PAOF::linearalgebra::{Vec2d};
use PAOF::commitment::{SisCommitment,KeyGen as CommKeyGen,Commitment,OpenComm};

fn main(){
    println!("Homomorphic commitment protocol based on SIS problem hardness.");
    let siscomm=SisCommitment::new(4099,50,50);
    let ck=siscomm.keygen();
    let m0=Vec2d{M:3,N:2,matrix:vec![vec![0,1,2],vec![3,4,5]]};
    let m1=Vec2d{M:3,N:2,matrix:vec![vec![0,1,2],vec![4,5,6]]};
    println!("M0: {:?}",m0);
    println!("M1: {:?}",m1);

    
    let (c0,r0)=siscomm.commitment(&ck,&m0);
    let (c1,r1)=siscomm.commitment(&ck,&m1);
    println!("C0: {:?}",c0);
    println!("C1: {:?}",c1);
    println!("Opening of C0+C1 given M0+M1:");
    let result=siscomm.open(&ck,&(&r0+&r1),&(&c0+&c1),&(&m0+&m1));
    println!("Result:{:?}",result);
}