# elina-rs

[![Crates.io](https://img.shields.io/crates/v/elina)](https://crates.io/crates/elina)

This library provides a safe Rust abstraction over the ELINA C library. To learn more about
ELINA, see their [GitHub](https://github.com/eth-sri/ELINA/) and [website](http://elina.ethz.ch/).

## Installation

You need ELINA installed to link this crate. Currently, you also need the ELINA header files
installed to build this crate, but this may change in the future.  
See [ELINA's website](http://elina.ethz.ch/) for instructions on how to install ELINA.

After you've installed ELINA, you can add this crate as a dependency in your `Cargo.toml`:
```toml
[dependencies]
elina = "0.1.1"
```

## Usage

Please read the [API documentation](https://docs.rs/elina).

See [`examples/showcase.rs`](examples/showcase.rs) for an example.  
You can run the example with `cargo run --example showcase`