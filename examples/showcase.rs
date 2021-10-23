use std::ops::Add;
use elina::ast::*;

fn main() {
    let env = Environment::new(vec!["x", "y", "z", "i"]);
    let man = OptPkManager::default();

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

    let meet_assn = meet.assign_copy(&man, &env, "z", &(y.clone() + Texpr::int(2)));
    println!("meet_assn (z = y + 2):");
    meet_assn.print(&man, &env);
    println!("meet_assn string: {}", meet_assn.to_string(&man, &env));
    println!("meet_assn satisfies z < 21: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(21))));
    println!("meet_assn satisfies z < 20: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(20))));
    println!("meet_assn satisfies z < 19: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(19))));
    println!("meet_assn satisfies z < 18: {}", meet_assn.satisfy(&man, &z.clone().lt(Texpr::int(18))));

    let x_bounds = meet_assn.get_bounds(&man, &env, "x");
    println!("x's bounds: {:?}", x_bounds);



    

    let x_gt_10 = Texpr::var(&env, "x").gt(Texpr::int(10));
    let x_lt_0 = Texpr::var(&env, "x").lt(Texpr::int(0));

    let top = Abstract::top(&man, &env);

    let hc_unsat = x_gt_10.clone().into_hcons().and(x_lt_0.clone().into());
    println!("Prev meet");
    let hc_unsat_meet = top.meet_copy(&man, &hc_unsat);
    println!("hc_unsat_meet:");
    hc_unsat_meet.print(&man, &env);

    let hc_or = hc_unsat.or(x.clone().lt(Texpr::int(1)).into());
    let hc_or_meet = top.meet_copy(&man, &hc_or);
    println!("hc_or_meet:");
    hc_or_meet.print(&man, &env);

    let hc_or_meet_joined = hc_or_meet.join_copy(&man, &meet_assn);
    println!("hc_or_meet_joined:");
    hc_or_meet_joined.print(&man, &env);

    let hc_or_meet_joined_string = hc_or_meet_joined.to_string(&man, &env);
    println!("hc_or_meet_joined string: {}", hc_or_meet_joined_string);

    // Widen stuff
    let x = x;
    let i = Texpr::var(&env, "i");
    let i_lt_x: Hcons = i.clone().lt(x.clone()).into();

    let mut top = Abstract::top(&man, &env);
    top.assign(&man, &env, "i", &y);
    top.meet(&man, &i_lt_x);

    let i0 = top.clone();
    top.assign(&man, &env, "i", &Texpr::int(1).add(y.clone()));
    top.meet(&man, &i_lt_x);
    let i1 = top.clone();
    println!("i0: {}", i0.to_string(&man, &env));
    println!("i1: {}", i1.to_string(&man, &env));
    println!("i0 widen i1: {}", (i0.widen_copy(&man, &i1)).to_string(&man, &env));


    // 'Testing' memory leaks

    // let mut meet_assn = meet_assn.clone();
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