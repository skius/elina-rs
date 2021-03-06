use std::ops::{Add, Mul};
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

    // meeting meet_assn with x < Texpr::top should change nothing
    let meet_assn2 = meet_assn.meet_copy(&man, &x.clone().lt(Texpr::top()));
    println!("meet_assn (with x < top): {}", meet_assn2.to_string(&man, &env));

    // neither should Hcons::top
    println!("meet_assn (meet with  top): {}", meet_assn.meet_copy(&man, &Hcons::Top).to_string(&man, &env));
    println!("meet_assn (meet with !top): {}", meet_assn.meet_copy(&man, &Hcons::Top.not()).to_string(&man, &env));

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


    let mut state = Abstract::top(&man, &env);
    state.assign(&man, &env, "y", &Texpr::int(2).mul(x.clone()));
    state.meet(&man, &(y.clone() + x.clone()).lt(Texpr::int(9)));
    println!("3x < 9: {}", state.to_string(&man, &env));
    println!("x bounds: {:?}", state.get_bounds(&man, &env, "x"));

    let mut state = Abstract::top(&man, &env);
    state.assign(&man, &env, "y", &Texpr::int(2).mul(x.clone()));
    println!("y = 2x: {}", state.to_string(&man, &env));
    state.add_dims(&man, 3, 1);
    state.assign_dim(&man, 4, &Texpr::int(5));
    println!("bounds of new state: {:?}", state.get_bounds_dim(&man, 4));
    println!("bounds of new state: {:?}", state.get_bounds_texpr(&man, &Texpr::int(100)));

    // Testing segfault
    // When one joins BOTTOM with something, where BOTTOM was obtained with unsat meet, segfault happens.
    // Segfault does not happen when Bottom is obtained with ::bottom()
    let mut state = hc_or_meet_joined.clone();
    let mut bot = state.clone();
    // let mut bot = Abstract::bottom(&man, &env);
    println!("bottom = state MEET false");
    bot.meet(&man, &Texpr::int(1).lt(Texpr::int(0)));
    println!("bottom JOIN state");
    bot.join(&man, &state);
    println!("{}", bot.to_string(&man, &env));
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