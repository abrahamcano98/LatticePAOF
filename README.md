# LatticePAOF
Polynomial arithmetic over finite fields for lattice-based encryption schemes

# Installing
### Via cargo
```
cargo install --git https://github.com/abrahamcano98/BTC-USD-ws
```
### Via git
```
git clone https://github.com/abrahamcano98/LatticePAOF
```


# Building
## Prerequisites
Rust and Cargo previously installed.

Then, execute (cargo only):
```
$ cargo build
```
# Basic usage
## Examples
The examples may be run by:
cargo run --example rlwe (Ring Learning With Errors encryption scheme)
cargo run --example commitment (Homomorphic commitment scheme)
cargo run --example zkpk (Amortized zero knowdledge proof of knowledge)
### Using traits, functions or structs in other files.
```use PAOF::mod::{element};```

For instance for using the function gauss_vec2d in the linearalgebra mod, type:

```use PAOF::linearalgebra::gauss_vec2d;```

