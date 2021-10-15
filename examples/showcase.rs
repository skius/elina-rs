use std::ops::Deref;
use elina_rs::ast::*;
use elina_sys::{ConsTyp, elina_interval_fprint, stdout};

fn main() {
    let env = Environment::new(vec!["x", "y", "z"]);
    let man = OptPkManager::default();

    let zero = Texpr::int(0);
    let x = Texpr::var(&env, "x");
    let y = Texpr::var(&env, "y");
    let z = Texpr::var(&env, "z");

    let upper = x.clone().lt(Texpr::int(20));
    let lower = x.clone().ge(Texpr::int(-10));

    let top = Abstract::top(&man, &env);
    println!("top:");
    top.print(&man, &env);
    println!("top satisfies x < 20: {}", top.satisfy(&man, &upper));
    println!();

    let mut meet = top.meet_copy(&man, &[&upper, &lower]);
    println!("meet:");
    meet.print(&man, &env);
    println!("meet satisfies x < 20: {}", meet.satisfy(&man, &upper));
    // println!("meet satisfies x2 < 400: {}", meet.satisfy(&man, &(x.clone() * x.clone()).lt(Texpr::int(400))));
    println!();

    let inter = y.clone().lt(x.clone());
    meet.meet(&man, &inter);
    println!("meet with `inter` constraint:");
    meet.print(&man, &env);
    println!("meet satisfies y < 20: {}", meet.satisfy(&man, &y.clone().lt(Texpr::int(20))));
    println!("meet satisfies y < 19: {}", meet.satisfy(&man, &y.clone().lt(Texpr::int(19))));
    println!("meet satisfies y < 18: {}", meet.satisfy(&man, &y.clone().lt(Texpr::int(18))));
    println!();

    let mut meet_assn = meet.assign_copy(&man, &env, "z", &(y.clone() + Texpr::int(2)));
    println!("meet_assn (z = y + 2):");
    meet_assn.print(&man, &env);
    println!("meet_assn satisfies z < 21: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(21))));
    println!("meet_assn satisfies z < 20: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(20))));
    println!("meet_assn satisfies z < 19: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(19))));
    println!("meet_assn satisfies z < 18: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(18))));

    let x_bounds = meet_assn.get_bounds(&man, &env, "x");
    println!("x's bounds: {:?}", x_bounds);


    // 'Testing' memory leaks

    // let mut i = 0;
    // loop {
    //     i += 1;
    //
    //     meet_assn.meet(&man, &[&upper, &lower]);
    //     if i % 10000 == 0 {
    //         println!("iter: {}", i);
    //         meet_assn.print(&man, &env);
    //     }
    // }
}