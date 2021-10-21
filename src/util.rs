// Could directly model these types in ast.rs and functions as methods on them,
// but for the time being this works

use std::ffi::CStr;
use std::os::raw::c_char;
use std::ptr::null_mut;

use elina_sys::{elina_lincons0_fprint, elina_lincons0_t, fclose, free, open_memstream};

use crate::ast::{Environment, Manager};

pub fn lincons0_to_string<M: Manager>(man: &M, env: &Environment, lincons0: *mut elina_lincons0_t) -> String {
    unsafe {
        let mut buf: *mut c_char = null_mut();
        let mut len = 0;
        let fd = open_memstream(&mut buf, &mut len);
        // to_env_names is okay here, because drop order guarantees it's only dropped after
        // the function has exited, at which point we don't need to strings anymore
        elina_lincons0_fprint(fd, lincons0, *env.to_env_names());
        fclose(fd);
        let res = CStr::from_ptr(buf).to_str().unwrap().to_string();
        free(buf as *mut _);
        return res;
        //
        //
        // UNFINISHED MANUAL IMPLEMENTATION:
        //
        // let linexpr = (*lincons0).linexpr0;
        // let mut res = linexpr0_to_string(man, env, linexpr);
        //
        // match (*lincons0).constyp.into() {
        //     ConsTyp::SUPEQ => res.push_str(" >= 0"),
        //     ConsTyp::SUP => res.push_str(" > 0"),
        //     ConsTyp::EQ => res.push_str(" = 0"),
        //     ConsTyp::DISEQ => res.push_str(" <> 0"),
        // }
        //
        // res
    }
}

// // https://adventures.michaelfbryan.com/posts/rust-closures-in-ffi/
// extern "C" fn each_test(i: size_t, dim: elina_dim_t, coeff: *mut elina_coeff_t) {
//     println!("foreach called with: {} {} {:?}", i, dim, coeff);
//
// }
//
// pub fn linexpr0_to_string<M: Manager>(man: &M, env: &Environment, linexpr0: *mut elina_linexpr0_t) -> String {
//     unsafe {
//         let fd = open_file(CString::new("temp.txt").unwrap().into_raw());
//         let env_names = env.to_env_names();
//         elina_linexpr0_fprint(fd, linexpr0, *env_names);
//         close_file(fd);
//
//         read_to_string("temp.txt").unwrap()
//         // foreach_linterm_of_linexpr0(linexpr0, Some(each_test))
//     }
// }

