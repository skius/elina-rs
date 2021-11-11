use std::ops::{Add, Mul};
use elina::ast::*;

fn main() {
    let env = Environment::new(vec!["x", "y", "z", "i"]);
    let man = OptPkManager::default();

    std::mem::drop(Abstract::bottom(&man, &env));
    println!("after ::bottom free");

    let mut top = Abstract::top(&man, &env);
    top.meet(&man, &Abstract::bottom(&man, &env));
    println!("after created bottom");
    std::mem::drop(top);
    println!("after created bottom free");
}