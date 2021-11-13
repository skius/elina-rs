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

### Projects

The following projects use `elina-rs` (feel free to open a PR with your own!):
- [skius/progge.rs](https://github.com/skius/progge.rs) - Program analysis playground for a simple, imperative language

## License

Licensed under either of

* Apache License, Version 2.0
  ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license
  ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.