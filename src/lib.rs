pub mod ast;

#[cfg(test)]
mod tests {
    use elina_sys::ConsTyp;
    use crate::ast::{Abstract, Environment, OptPkManager, Tcons, Texpr};

    #[test]
    fn integration_test() {
        let env = Environment::new(vec!["x"]);
        let leaf_5 = Texpr::int(5);
        let leaf_x = Texpr::var(&env, "x");
        let x_plus_5 = leaf_x + leaf_5;

        let tcons1 = Tcons::new(x_plus_5, ConsTyp::SUPEQ);

        let tcons2 = Tcons::new(Texpr::int(20) - Texpr::var(&env, "x"), ConsTyp::SUPEQ);

        let man = OptPkManager::default();

        let top = Abstract::top(&man, &env);
        let meet_res = top.meet_copy(&man, &tcons1);
        let meet_res = meet_res.meet_copy(&man, &meet_res);
        let meet_res = meet_res.meet_copy(&man, &tcons2);

        meet_res.print(&man, &env);
    }

    #[test]
    fn new_environment() {
        let env = Environment::new(vec!["hello"]);
        let leaf = Texpr::var(&env, "hello");

        println!("{:?}", leaf)
    }


    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
