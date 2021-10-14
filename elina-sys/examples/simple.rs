use std::ffi::CString;
use std::thread::sleep;
use std::time::Duration;
use elina_sys::*;

fn main() {
    unsafe {
        let intdim = 5;
        let realdim = 0;
        let mut names_sl = ["x0", "x1", "x2", "x3", "x4"].map(|s| {
            let mut cstring = CString::new(s);
            cstring.unwrap().into_raw()
        });
        let names = names_sl.as_mut_ptr();
        std::mem::forget(names_sl);
        // println!("{:?}", CString::from_raw(*names));

        let man = opt_pk_manager_alloc(false_);

        let x0 = elina_texpr0_dim(0);
        let cst_5 = elina_texpr0_cst_scalar_int(5);
        let x0_plus_5 = elina_texpr0_binop(elina_texpr_op_t_ELINA_TEXPR_ADD, x0, cst_5, elina_texpr_rtype_t_ELINA_RTYPE_INT, elina_texpr_rdir_t_ELINA_RDIR_ZERO);
        let top = elina_abstract0_top(man, intdim, realdim);


        println!("before print texpr0");
        elina_texpr0_print(x0_plus_5, names);
        let _ = *x0_plus_5;
        println!("after print texpr0");

        let mut tcons0 = elina_tcons0_t {
            texpr0: x0_plus_5,
            constyp: elina_constyp_t_ELINA_CONS_SUPEQ,
            scalar: 0 as *mut elina_scalar_t,
        };

        elina_tcons0_fprint(stdout, &mut tcons0, names);

        let mut tconsarr = elina_tcons0_array_make(1);
        *tconsarr.p = tcons0;


        let meet_with_tcons = elina_abstract0_meet_tcons_array(man , false_, top, &mut tconsarr);

        let mut lcarr = elina_abstract0_to_lincons_array(man, meet_with_tcons);
        elina_lincons0_array_print(&mut lcarr, names);
        println!("after tcons meet stuff");

        let res = elina_abstract0_sat_tcons(man, meet_with_tcons, &mut tcons0);
        println!("{}", if res == true_ {true} else {false});

        // so far so good
        return;

        // opt_pk_meet_tcons_array(man, false_, , &mut tconsarr);

        let linexpr = elina_linexpr0_alloc(elina_linexpr_discr_t_ELINA_LINEXPR_SPARSE, 0);
        elina_scalar_set_to_int((*linexpr).cst.val.scalar, 5, elina_scalar_discr_t_ELINA_SCALAR_MPQ);

        println!("before other");
        // let linexpr = elina_intlinearize_texpr0_intlinear(man, x0_plus_5, elina_scalar_discr_t_ELINA_SCALAR_MPQ);
        // //this line crashes:
        // dbg!((*linexpr).size);
        //
        // // linexpr0_all
        //
        // elina_scalar_set_to_int((*linexpr).cst.val.scalar, 0, elina_scalar_discr_t_ELINA_SCALAR_MPQ);
        println!("before print linexpr");
        elina_linexpr0_print(linexpr, names);
        println!("after print linexpr");


        let scalar = elina_scalar_alloc();
        elina_scalar_set_int(scalar, 0);
        //
        let mut cons = elina_lincons0_make(elina_constyp_t_ELINA_CONS_SUPEQ, linexpr, scalar);
        println!("Test");

        elina_lincons0_print(&mut cons, names);
        println!("Test");
        let mut consarr = elina_lincons0_array_make(1);
        *consarr.p = cons;

        let meet_top = elina_abstract0_meet_lincons_array(man, false_, top, &mut consarr);



        println!("first thing: ");
        let mut lcarr = elina_abstract0_to_lincons_array(man, meet_top);
        elina_lincons0_array_print(&mut lcarr, names);

        // elina_abstract0_assign_texpr()
        //
        // let tcons = elina_tcons0_t {}

        // let meet = elina_abstract0_meet(man, false_, top, cons);
        let names = ["x0", "x1", "x2", "x3", "x4"].map(|s| {
            let mut cstring = CString::new(s);
            cstring.unwrap().into_raw()
        }).as_mut_ptr();
        println!("names: {:?}", names);
        elina_abstract0_fprint(stdout, man, top, names);
        println!("after: ");
        let mut lcarr = elina_abstract0_to_lincons_array(man, top);
        elina_lincons0_array_print(&mut lcarr, names);

        //sleep(Duration::from_millis(1000));

    }
}