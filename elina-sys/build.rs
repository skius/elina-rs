fn main() {
    if let Ok(_) = std::env::var("DOCS_RS") {
        // We don't need to link for building the docs.
        return;
    }


    // Tell cargo to tell rustc to link the system ELINA
    // shared library.
    println!("cargo:rustc-link-lib=elinalinearize");
    println!("cargo:rustc-link-lib=elinaux");
    println!("cargo:rustc-link-lib=optpoly");
    println!("cargo:rustc-link-lib=gmp");
    println!("cargo:rustc-link-lib=mpfr");

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-changed=wrapper.c");



    cc::Build::new()
        .file("wrapper.c")
        .compile("wrapper");

    // // The bindgen::Builder is the main entry point
    // // to bindgen, and lets you build up options for
    // // the resulting bindings.
    // let bindings = bindgen::Builder::default()
    //     // The input header we would like to generate
    //     // bindings for.
    //     .header("wrapper.h")
    //     // Tell cargo to invalidate the built crate whenever any of the
    //     // included header files changed.
    //     .parse_callbacks(Box::new(bindgen::CargoCallbacks))
    //
    //     // Blocklist the specific ELINA consts that cause issues
    //     .blocklist_item("FP_NAN")
    //     .blocklist_item("FP_INFINITE")
    //     .blocklist_item("FP_ZERO")
    //     .blocklist_item("FP_SUBNORMAL")
    //     .blocklist_item("FP_NORMAL")
    //
    //     // Finish the builder and generate the bindings.
    //     .generate()
    //     // Unwrap the Result and panic on failure.
    //     .expect("Unable to generate bindings");
    //
    // // Write the bindings to the $OUT_DIR/bindings.rs file.
    // let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    // bindings
    //     .write_to_file(out_path.join("bindings.rs"))
    //     .expect("Couldn't write bindings!");
}