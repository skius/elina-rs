#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[repr(u32)]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum TexprUnop {
    Neg = elina_texpr_op_t_ELINA_TEXPR_NEG,
}

#[repr(u32)]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum TexprBinop {
    Add = elina_texpr_op_t_ELINA_TEXPR_ADD,
    Sub = elina_texpr_op_t_ELINA_TEXPR_SUB,
    Mul = elina_texpr_op_t_ELINA_TEXPR_MUL,
    Div = elina_texpr_op_t_ELINA_TEXPR_DIV,
    Mod = elina_texpr_op_t_ELINA_TEXPR_MOD,
}

#[repr(u32)]
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum ConsTyp {
    SUPEQ = elina_constyp_t_ELINA_CONS_SUPEQ,
    SUP = elina_constyp_t_ELINA_CONS_SUP,
    EQ = elina_constyp_t_ELINA_CONS_EQ,
    DISEQ = elina_constyp_t_ELINA_CONS_DISEQ,
    //MODEQ, needs support in Tcons in the form of a proper scalar
}

pub fn bool_from_c_bool(b: bool_) -> bool {
    b == true_
}

pub fn c_bool_from_bool(b: bool) -> bool_ {
    if b {
        true_
    } else {
        false_
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;
    use crate::{elina_abstract0_approximate, elina_abstract0_fprint, elina_abstract0_meet, elina_abstract0_top, elina_constyp_t_ELINA_CONS_SUPEQ, elina_manager_alloc, elina_tcons0_t, elina_texpr0_cst_scalar_int, elina_texpr0_dim, elina_texpr0_print, false_, opt_pk_manager_alloc, stderr, stdout};



    #[test]
    fn manager_alloc() {
        unsafe {
            let intdim = 5;
            let realdim = 0;


            let man = opt_pk_manager_alloc(false_);

            let x0 = elina_texpr0_dim(0);
            let cst_5 = elina_texpr0_cst_scalar_int(5);
            let top = elina_abstract0_top(man, intdim, realdim);

            // let meet = elina_abstract0_meet(man, false_, top, cons);
            let names = ["x0", "x1", "x2", "x3", "x4"].map(|s| {
                let mut cstring = CString::new(s);
                cstring.unwrap().into_raw()
            }).as_mut_ptr();
            let ret = elina_abstract0_fprint(stderr, man, top, names);


            // elina_texpr0_print()

        }
    }

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
