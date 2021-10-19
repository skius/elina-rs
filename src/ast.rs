use std::borrow::Borrow;
use std::collections::HashMap;
use std::ffi::{CStr, CString};
use std::fmt::{Debug, Formatter};
use std::ptr::null_mut;

pub use elina_sys::{ConsTyp, TexprBinop, TexprUnop};
use elina_sys::{__gmpq_get_str, __gmpz_export, bool_from_c_bool, c_bool_from_bool, elina_abstract0_assign_texpr, elina_abstract0_bottom, elina_abstract0_bound_dimension, elina_abstract0_copy, elina_abstract0_free, elina_abstract0_is_bottom, elina_abstract0_is_top, elina_abstract0_join, elina_abstract0_meet, elina_abstract0_meet_tcons_array, elina_abstract0_sat_tcons, elina_abstract0_t, elina_abstract0_to_lincons_array, elina_abstract0_top, elina_constyp_t, elina_constyp_t_ELINA_CONS_DISEQ, elina_constyp_t_ELINA_CONS_EQ, elina_constyp_t_ELINA_CONS_SUPEQ, elina_dim_t, elina_interval_free, elina_lincons0_array_clear, elina_lincons0_array_print, elina_manager_free, elina_manager_t, elina_scalar_free, elina_scalar_t, elina_tcons0_array_make, elina_tcons0_t, elina_texpr0_binop, elina_texpr0_copy, elina_texpr0_cst_scalar_int, elina_texpr0_dim, elina_texpr0_free, elina_texpr0_t, elina_texpr0_unop, elina_texpr_op_t, elina_texpr_rdir_t_ELINA_RDIR_ZERO, elina_texpr_rtype_t_ELINA_RTYPE_INT, false_, free, opt_pk_manager_alloc, true_};

/// Provides the implementations of different abstract domains.
pub trait Manager {
    /// Returns the `elina_manager_t` pointer required for internal function calls.
    unsafe fn as_manager_ptr(&self) -> *mut elina_manager_t;
}

/// ELINA's polyhedra domain manager.
///
/// Wraps `elina_manager_t` obtained by `opt_pkt_manager_alloc`.
pub struct OptPkManager {
    elina_manager: *mut elina_manager_t,
}

impl Default for OptPkManager {
    fn default() -> Self {
        unsafe {
            OptPkManager {
                elina_manager: opt_pk_manager_alloc(false_),
            }
        }
    }
}

impl Drop for OptPkManager {
    fn drop(&mut self) {
        unsafe {
            elina_manager_free(self.elina_manager);
        }
    }
}

impl Manager for OptPkManager {
    unsafe fn as_manager_ptr(&self) -> *mut elina_manager_t {
        self.elina_manager
    }
}

/// Keeps track of variables names and their associated dimension in the abstract element.
#[derive(Debug)]
pub struct Environment {
    var_to_dim: HashMap<String, elina_dim_t>,
}

impl Environment {
    /// Returns an `Environment` consisting of the variables in the ordered vector `int_names`.
    pub fn new<S>(int_names: Vec<S>) -> Environment
    where
        S: AsRef<str>,
    {
        Environment {
            var_to_dim: int_names
                .into_iter()
                .enumerate()
                .map(|(i, v)| (v.as_ref().to_owned(), i as u32))
                .collect(),
        }
    }
}

// pub trait Texpr1 {
//     fn as_texpr0_ptr(&mut self, env: &Environment) -> *mut elina_texpr0_t;
// }
//
// pub enum Texpr11 {
//     Unop(),
//     Binop(TexprBinop, Box<>),
// }

/// A tree-based expression.
///
/// Stores its expression in the usual tree-based format from e.g. abstract syntax trees. This
/// allows for simple clients of this type.
///
/// Wraps `elina_texpr0_t`.
#[derive(Debug)]
pub struct Texpr {
    elina_texpr0: *mut elina_texpr0_t,
}

impl Texpr {
    /// Returns a `Texpr` representing a constant integer.
    pub fn int(i: i64) -> Texpr {
        unsafe {
            Texpr {
                elina_texpr0: elina_texpr0_cst_scalar_int(i),
            }
        }
    }

    /// Returns a `Texpr` representing a variable with name `s`.
    ///
    /// Requires `env`, because internally variables are simply dimensions of the abstract element.
    pub fn var<S: Borrow<str>>(env: &Environment, s: S) -> Texpr {
        unsafe {
            Texpr {
                elina_texpr0: elina_texpr0_dim(env.var_to_dim[s.borrow()]),
            }
        }
    }

    /// Returns a `Texpr` representing the variable associated with dimension `dim`.
    pub fn var_dim(dim: u32) -> Texpr {
        unsafe {
            Texpr {
                elina_texpr0: elina_texpr0_dim(dim),
            }
        }
    }

    /// Returns a `Texpr` representing the binary expression `left op right`.
    pub fn binop(op: TexprBinop, left: Texpr, right: Texpr) -> Texpr {
        unsafe {
            let res = Texpr {
                elina_texpr0: elina_texpr0_binop(
                    op as elina_texpr_op_t,
                    left.elina_texpr0,
                    right.elina_texpr0,
                    elina_texpr_rtype_t_ELINA_RTYPE_INT,
                    elina_texpr_rdir_t_ELINA_RDIR_ZERO,
                ),
            };

            std::mem::forget(left);
            std::mem::forget(right);

            res
        }
    }

    /// Returns a `Texpr` representing the unary expression `op inner`.
    pub fn unop(op: TexprUnop, inner: Texpr) -> Texpr {
        unsafe {
            let res = Texpr {
                elina_texpr0: elina_texpr0_unop(
                    op as elina_texpr_op_t,
                    inner.elina_texpr0,
                    elina_texpr_rtype_t_ELINA_RTYPE_INT,
                    elina_texpr_rdir_t_ELINA_RDIR_ZERO,
                ),
            };

            std::mem::forget(inner);

            res
        }
    }

    /*
        Tcons constraint generating operators
        These have a special meaning in Rust and cannot be overloaded
    */

    pub fn lt(self, rhs: Texpr) -> Tcons {
        /*
            self < rhs
            0 < rhs - self
            0 <= rhs - self - 1
        */
        Tcons::new(rhs - self - Texpr::int(1), ConsTyp::SUPEQ)
    }

    pub fn le(self, rhs: Texpr) -> Tcons {
        /*
            self <= rhs
            0 <= rhs - self
        */
        Tcons::new(rhs - self, ConsTyp::SUPEQ)
    }

    pub fn gt(self, rhs: Texpr) -> Tcons {
        /*
            self > rhs
            self - rhs > 0
            self - rhs - 1 >= 0
        */
        Tcons::new(self - rhs - Texpr::int(1), ConsTyp::SUPEQ)
    }

    pub fn ge(self, rhs: Texpr) -> Tcons {
        /*
            self >= rhs
            self - rhs >= 0
        */
        Tcons::new(self - rhs, ConsTyp::SUPEQ)
    }

    pub fn _eq(self, rhs: Texpr) -> Tcons {
        /*
            self = rhs
            self - rhs = 0
        */
        Tcons::new(self - rhs, ConsTyp::EQ)
    }

    pub fn _ne(self, rhs: Texpr) -> Tcons {
        /*
            self != rhs
            self - rhs != 0
        */
        // TODO: Possible that it's better to treat this as join of < and >
        Tcons::new(self - rhs, ConsTyp::DISEQ)
    }
}

impl std::ops::Add<Texpr> for Texpr {
    type Output = Texpr;

    fn add(self, rhs: Texpr) -> Self::Output {
        Texpr::binop(TexprBinop::Add, self, rhs)
    }
}

impl std::ops::Sub<Texpr> for Texpr {
    type Output = Texpr;

    fn sub(self, rhs: Texpr) -> Self::Output {
        Texpr::binop(TexprBinop::Sub, self, rhs)
    }
}

impl std::ops::Mul<Texpr> for Texpr {
    type Output = Texpr;

    fn mul(self, rhs: Texpr) -> Self::Output {
        Texpr::binop(TexprBinop::Mul, self, rhs)
    }
}

impl std::ops::Div<Texpr> for Texpr {
    type Output = Texpr;

    fn div(self, rhs: Texpr) -> Self::Output {
        Texpr::binop(TexprBinop::Div, self, rhs)
    }
}

impl std::ops::Rem<Texpr> for Texpr {
    type Output = Texpr;

    fn rem(self, rhs: Texpr) -> Self::Output {
        Texpr::binop(TexprBinop::Mod, self, rhs)
    }
}

impl std::ops::Neg for Texpr {
    type Output = Texpr;

    fn neg(self) -> Self::Output {
        Texpr::unop(TexprUnop::Neg, self)
    }
}

impl Clone for Texpr {
    fn clone(&self) -> Self {
        unsafe {
            Texpr {
                elina_texpr0: elina_texpr0_copy(self.elina_texpr0),
            }
        }
    }
}

impl Drop for Texpr {
    fn drop(&mut self) {
        unsafe {
            elina_texpr0_free(self.elina_texpr0);
        }
    }
}

#[derive(Debug, Clone)]
pub enum HconsBinop {
    Or,
    And,
}

impl HconsBinop {
    pub fn negation(&self) -> HconsBinop {
        use HconsBinop::*;

        match self {
            Or => And,
            And => Or,
        }
    }
}

#[derive(Debug, Clone)]
pub enum HconsUnop {
    Not,
}

/// A multi-level constraint.
#[derive(Clone)]
pub enum Hcons {
    /// Wraps a [`Tcons`].
    Leaf(Tcons),
    /// A binary operation on two constraints.
    Binop(HconsBinop, Box<Hcons>, Box<Hcons>),
    /// A unary operation on a constraint.
    Unop(HconsUnop, Box<Hcons>),
}

impl Hcons {
    pub fn and(left: Hcons, right: Hcons) -> Hcons {
        Hcons::Binop(
            HconsBinop::And,
            Box::new(left),
            Box::new(right)
        )
    }

    pub fn or(left: Hcons, right: Hcons) -> Hcons {
        Hcons::Binop(
            HconsBinop::Or,
            Box::new(left),
            Box::new(right)
        )
    }

    pub fn not(inner: Hcons) -> Hcons {
        Hcons::Unop(HconsUnop::Not, Box::new(inner))
    }

    pub fn leaf(tcons: Tcons) -> Hcons {
        Hcons::Leaf(tcons)
    }

    pub fn negation(&self) -> Hcons {
        use Hcons::*;

        match self {
            Leaf(tcons) => Leaf(tcons.negation()),
            Binop(bop, left, right) =>
                Binop(
                    bop.negation(),
                    Box::new(left.negation()),
                    Box::new(right.negation())
                ),
            Unop(HconsUnop::Not, inner) => *inner.clone(),
        }
    }
}

/// A tree-based constraint.
///
/// Wraps `elina_tcons0_t`.
pub struct Tcons {
    elina_tcons0: *mut elina_tcons0_t,
}

impl Tcons {
    /// Creates a new `Tcons`, consuming `texpr` in the process.
    ///
    /// Alternatively, you can use the various [`Texpr::lt`], [`Texpr::gt`], etc. functions.
    pub fn new(texpr: Texpr, cons_typ: ConsTyp) -> Tcons {
        unsafe {
            let res = Tcons::from_raw(texpr.elina_texpr0, cons_typ as elina_constyp_t);

            std::mem::forget(texpr);

            res
        }
    }

    unsafe fn from_raw(texpr0: *mut elina_texpr0_t, cons_typ: elina_constyp_t) -> Tcons {
        let inner_tcons = Box::new(elina_tcons0_t {
            texpr0,
            constyp: cons_typ,
            scalar: 0 as *mut elina_scalar_t, // doesn't support EQMOD
        });
        let res = Tcons {
            elina_tcons0: Box::into_raw(inner_tcons),
        };
        res
    }

    pub fn negation(&self) -> Tcons {
        unsafe {
            let cloned = self.clone();
            match (*self.elina_tcons0).constyp {
                x if x == elina_constyp_t_ELINA_CONS_EQ => {
                    (*cloned.elina_tcons0).constyp = ConsTyp::DISEQ as elina_constyp_t;
                    cloned
                },
                x if x == elina_constyp_t_ELINA_CONS_DISEQ => {
                    (*cloned.elina_tcons0).constyp = ConsTyp::EQ as elina_constyp_t;
                    cloned
                }
                x if x == elina_constyp_t_ELINA_CONS_SUPEQ => {
                    let texpr = Texpr { elina_texpr0: (*cloned.elina_tcons0).texpr0 };
                    // std::mem::forget(cloned);

                    // not (texpr >= 0)
                    // texpr < 0
                    // -texpr > 0
                    // -texpr - 1 >= 0

                    let new_texpr = Texpr::int(-1) * texpr - Texpr::int(1);
                    (*cloned.elina_tcons0).texpr0 = new_texpr.elina_texpr0;
                    std::mem::forget(new_texpr);
                    cloned
                }
                _ => todo!()
            }
        }
    }
}

impl Clone for Tcons {
    fn clone(&self) -> Self {
        unsafe {
            Tcons::from_raw(
                elina_texpr0_copy((*self.elina_tcons0).texpr0),
                (*self.elina_tcons0).constyp,
            )
        }
    }
}

impl Drop for Tcons {
    // TODO: idea, add NULL pointer check to all Drop, then we don't need to use as many mem::forget's
    fn drop(&mut self) {
        unsafe {
            let tcons0 = Box::from_raw(self.elina_tcons0);

            elina_texpr0_free(tcons0.texpr0);
            if tcons0.scalar as u64 != 0 {
                elina_scalar_free(tcons0.scalar);
            }
        }
    }
}

///  An element of the abstract domain lattice.
///
/// In ELINA, a single `Abstract` models mappings from variables to sets of numbers.
///
/// Wraps `elina_abstract0_t`.
pub struct Abstract {
    elina_abstract0: *mut elina_abstract0_t,
}

impl Abstract {
    /// Returns a new `Abstract` element representing Top (⊤) in the lattice.
    pub fn top<M: Manager>(man: &M, env: &Environment) -> Abstract {
        unsafe {
            Abstract {
                elina_abstract0: elina_abstract0_top(
                    man.as_manager_ptr(),
                    env.var_to_dim.len() as u64,
                    0,
                ),
            }
        }
    }

    /// Returns a new `Abstract` element representing Bottom (⊥) in the lattice.
    pub fn bottom<M: Manager>(man: &M, env: &Environment) -> Abstract {
        unsafe {
            Abstract {
                elina_abstract0: elina_abstract0_bottom(
                    man.as_manager_ptr(),
                    env.var_to_dim.len() as u64,
                    0,
                ),
            }
        }
    }

    /// Returns `true` if `self` satisfies `tcons`, i.e. `self` ⊆ `tcons`.
    pub fn satisfy<M: Manager>(&self, man: &M, tcons: &Tcons) -> bool {
        unsafe {
            bool_from_c_bool(elina_abstract0_sat_tcons(
                man.as_manager_ptr(),
                self.elina_abstract0,
                tcons.elina_tcons0,
            ))
        }
    }

    /// Performs the meet operation on the lattice with `self` and `other`, and stores the
    /// result in `self`.
    ///
    /// See the copying counterpart at [`Abstract::meet_copy`].
    pub fn meet<M: Manager, MT: Meetable + ?Sized>(&mut self, man: &M, other: &MT) {
        other.meet_with(man, self);
    }

    /// Returns the result of the meet operation on the lattice with `self` and `other`.
    ///
    /// See the mutating counterpart at [`Abstract::meet`].
    pub fn meet_copy<M: Manager, MT: Meetable + ?Sized>(&self, man: &M, other: &MT) -> Abstract {
        other.meet_with_copy(man, self)
    }

    /// Performs the join operation on the lattice with `self` and `other`, and stores the
    /// result in `self`.
    ///
    /// See the copying counterpart at [`Abstract::join_copy`].
    pub fn join<M: Manager, JT: Joinable + ?Sized>(&mut self, man: &M, other: &JT) {
        other.join_with(man, self);
    }

    /// Returns the result of the join operation on the lattice with `self` and `other`.
    ///
    /// See the mutating counterpart at [`Abstract::join`].
    pub fn join_copy<M: Manager, JT: Joinable + ?Sized>(&self, man: &M, other: &JT) -> Abstract {
        other.join_with_copy(man, self)
    }

    /// Assigns `var` to `texpr` in `self`.
    ///
    /// This function can be used to model mutable variables.
    ///
    /// See the copying counterpart at [`Abstract::assign_copy`].
    pub fn assign<M: Manager, S: Borrow<str>>(
        &mut self,
        man: &M,
        env: &Environment,
        var: S,
        texpr: &Texpr,
    ) {
        unsafe {
            elina_abstract0_assign_texpr(
                man.as_manager_ptr(),
                true_,
                self.elina_abstract0,
                env.var_to_dim[var.borrow()],
                texpr.elina_texpr0,
                std::ptr::null_mut(),
            );
        }
    }

    /// Returns a new `Abstract` representing `self` after `var` has been assigned `texpr`.
    ///
    /// This function can be used to model mutable variables.
    ///
    /// See the mutating counterpart at [`Abstract::assign`].
    pub fn assign_copy<M, S>(&self, man: &M, env: &Environment, var: S, texpr: &Texpr) -> Abstract
    where
        M: Manager,
        S: Borrow<str>,
    {
        unsafe {
            Abstract {
                elina_abstract0: elina_abstract0_assign_texpr(
                    man.as_manager_ptr(),
                    false_,
                    self.elina_abstract0,
                    env.var_to_dim[var.borrow()],
                    texpr.elina_texpr0,
                    std::ptr::null_mut(),
                ),
            }
        }
    }

    /// Returns `true` if `self` is Top.
    pub fn is_top<M: Manager>(&self, man: &M) -> bool {
        unsafe {
            bool_from_c_bool(elina_abstract0_is_top(
                man.as_manager_ptr(),
                self.elina_abstract0
            ))
        }
    }

    /// Returns `true` if `self` is Bottom.
    pub fn is_bottom<M: Manager>(&self, man: &M) -> bool {
        unsafe {
            bool_from_c_bool(elina_abstract0_is_bottom(
                man.as_manager_ptr(),
                self.elina_abstract0
            ))
        }
    }

    /// Returns the bounds of variable `var` in `self`.
    pub fn get_bounds<M, S>(&self, man: &M, env: &Environment, var: S) -> Interval
    where
        M: Manager,
        S: Borrow<str>,
    {
        unsafe {
            // the OPT_PK manager does not support getting a texpr's bound, thus we need a workaround
            // let interval_ptr = elina_abstract0_bound_texpr(man.as_manager_ptr(), self.elina_abstract0, texpr.elina_texpr0);

            // let linexpr = elina_intlinearize_texpr0_intlinear(
            //     man.as_manager_ptr(),
            //     texpr.elina_texpr0,
            //     elina_scalar_discr_t_ELINA_SCALAR_MPQ
            // );
            // let lin = *linexpr;
            // println!("{:?}", &lin as *const elina_linexpr0_t);
            //
            // let interval_ptr = elina_abstract0_bound_linexpr(
            //     man.as_manager_ptr(),
            //     self.elina_abstract0,
            //     linexpr
            // );

            let dim = env.var_to_dim[var.borrow()];
            let interval_ptr =
                elina_abstract0_bound_dimension(man.as_manager_ptr(), self.elina_abstract0, dim);

            let inf = *(*interval_ptr).inf;
            let inf = *inf.val.mpq;
            let inf_denom = inf._mp_den;
            let inf_enum = inf._mp_num;

            let mut result_denom = 0u64;
            let mut result_enum = 0u64;
            __gmpz_export(
                &mut result_denom as *mut _ as _,
                null_mut(),
                -1,
                8,
                0,
                0,
                &inf_denom,
            );
            __gmpz_export(
                &mut result_enum as *mut _ as _,
                null_mut(),
                -1,
                8,
                0,
                0,
                &inf_enum,
            );

            let inf_str_ptr = __gmpq_get_str(null_mut(), 10, &inf);
            let inf_str = CStr::from_ptr(inf_str_ptr);
            // println!("inf_str: {:?}", inf_str);

            let inf = if result_denom == 0 {
                Bound::NegInfinity
            } else {
                Bound::Num(inf_str.to_str().unwrap().parse().unwrap())
            };
            free(inf_str_ptr as *mut _);

            // println!("inf: {:?}", inf);

            let sup = *(*interval_ptr).sup;
            let sup = *sup.val.mpq;
            let sup_denom = sup._mp_den;
            let sup_enum = sup._mp_num;

            let mut result_denom = 0u64;
            let mut result_enum = 0u64;
            __gmpz_export(
                &mut result_denom as *mut _ as _,
                null_mut(),
                -1,
                8,
                0,
                0,
                &sup_denom,
            );
            __gmpz_export(
                &mut result_enum as *mut _ as _,
                null_mut(),
                -1,
                8,
                0,
                0,
                &sup_enum,
            );

            let sup_str_ptr = __gmpq_get_str(null_mut(), 10, &sup);
            let sup_str = CStr::from_ptr(sup_str_ptr);
            // println!("sup_str: {:?}", sup_str);

            let sup = if result_denom == 0 {
                Bound::PosInfinity
            } else {
                Bound::Num(sup_str.to_str().unwrap().parse().unwrap())
            };
            free(sup_str_ptr as *mut _);

            // println!("sup: {:?}", sup);

            // gives f64
            // let inf = __gmpq_get_d((*inf).val.mpq);
            // println!("{}", inf);

            elina_interval_free(interval_ptr);

            Interval(inf, sup)
        }
    }

    /// Prints `self` to `stdout`.
    pub fn print<M: Manager>(&self, man: &M, env: &Environment) {
        unsafe {
            let rev_env = env
                .var_to_dim
                .iter()
                .map(|(k, v)| (v.to_owned(), k.to_owned()))
                .collect::<HashMap<elina_dim_t, String>>();

            let mut names = (0..env.var_to_dim.len() as u32)
                .into_iter()
                .map(|d| CString::new(rev_env[&d].as_str()).unwrap().into_raw())
                .collect::<Vec<_>>();

            let names_ptr = names.as_mut_ptr();
            // std::mem::forget(names);

            let mut lincons_arr =
                elina_abstract0_to_lincons_array(man.as_manager_ptr(), self.elina_abstract0);
            // println!("Reached");
            // println!("{:?}", CString::from_raw(*names_ptr));
            // println!("{:?}", lincons_arr.p);
            elina_lincons0_array_print(&mut lincons_arr, names_ptr);

            elina_lincons0_array_clear(&mut lincons_arr);
            std::mem::drop(lincons_arr);

            names
                .into_iter()
                .for_each(|ptr| std::mem::drop(CString::from_raw(ptr)));
        }
    }
}

impl Drop for Abstract {
    fn drop(&mut self) {
        unsafe {
            if bool_from_c_bool(elina_abstract0_is_bottom((*self.elina_abstract0).man, self.elina_abstract0)) {
                // println!("Abstract drop ignored because it's bottom");
                return;
            }
            // println!("Drop incoming:");
            // let man_ptr = (*self.elina_abstract0).man;
            // println!("Dropping Abstract:!!!");
            // let man = *man_ptr;
            // println!("{:?}", &man as *const _);
            elina_abstract0_free((*self.elina_abstract0).man, self.elina_abstract0);
            // println!("dropped.");
        }
    }
}

impl Clone for Abstract {
    fn clone(&self) -> Self {
        unsafe {
            Abstract {
                elina_abstract0: elina_abstract0_copy(
                    (*self.elina_abstract0).man,
                    self.elina_abstract0,
                ),
            }
        }
    }
}

// impl Display for Abstract {
//     fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
//         unsafe {
//             let man = self.elina_abstract0.man;
//             let lincons_arr = elina_abstract0_to_lincons_array(man, self.elina_abstract0);
//
//             let lincons_sl = slice_from_raw_parts(lincons_arr.p, lincons_arr.size as usize);
//
//             f.write_str("{ ");
//             for lincons in &*lincons_sl {
//                 f.write_str(lincons.linexpr0.discr)
//             }
//
//             Ok(())
//         }
//     }
// }

// TODO: Add docs, examples
// also add a custom struct that allows AND/ORs between Tcons by lazily storing them
// and only evaluating them when Meet/Join is called

/// An element that [`Abstract`] can meet with.
///
/// Outside this trait, member functions are only called by
/// [`Abstract::meet`] and [`Abstract::meet_copy`].
///
/// # Implementation
///
/// Implementors of this trait must either implement [`Meetable::meet_internal`] or both
/// [`Meetable::meet_with`] and [`Meetable::meet_with_copy`].
pub trait Meetable {
    /// Returns a pointer to the internal meet result of `self` with `other`.
    ///
    /// If `destructive` is `true`, `other` is internally mutated and its pointer is returned,
    /// if `destructive` is `false`, `other` is not mutated and a pointer to the newly allocated
    /// internal `elina_abstract0_t` is returned.
    ///
    /// # Panics
    ///
    /// This function is only called by [`Meetable::meet_with`] and [`Meetable::meet_with_copy`].
    /// It panics when neither `meet_internal` nor the two `meet_with`'s are implemented.
    ///
    /// # Safety
    ///
    /// This is unsafe, because `other` is internally mutated exactly if `destructive` is true.
    /// The caller must ensure that they have mutable permissions for `other` when calling this
    /// function with `destructive` set to true.
    #[allow(unused_variables)]
    unsafe fn meet_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        panic!("Meetable::meet_internal's definition must be overridden if you don't provide implementations for meet_with and meet_with_copy");
    }

    // reverse, self.meet_with(other) = elina_abstract0_meet...(other, self)
    fn meet_with<M: Manager>(&self, man: &M, other: &mut Abstract) {
        unsafe {
            // Mutates "other"
            other.elina_abstract0 = self.meet_internal(man, other.elina_abstract0, true);
        }
    }
    fn meet_with_copy<M: Manager>(&self, man: &M, other: &Abstract) -> Abstract {
        unsafe {
            let new_abs_ptr = self.meet_internal(man, other.elina_abstract0, false);
            Abstract {
                elina_abstract0: new_abs_ptr,
            }
        }
    }
}

impl Meetable for Hcons {
    unsafe fn meet_internal<M: Manager>(&self, man: &M, other: *mut elina_abstract0_t, destructive: bool) -> *mut elina_abstract0_t {
        use Hcons::*;
        use HconsBinop::*;
        use HconsUnop::*;

        match self {
            Leaf(tcons) => tcons.meet_internal(man, other, destructive),
            Binop(And, left, right) => {
                if destructive {
                    let res = left.meet_internal(man, other, true);
                    let res = right.meet_internal(man, res, true);
                    res
                } else {
                    let interim = left.meet_internal(man, other, false);
                    let res = right.meet_internal(man, interim, true);
                    res
                }
                // // first do `destructive` meet
                // let left_res = left.meet_internal(man, other, destructive);
                // let mut left_abs = Abstract { elina_abstract0: left_res };
                // // then do mutating meet on result
                // println!("reached binop");
                // dbg!(destructive);
                // // let meet_ptr = right.meet_internal(man, &left_abs, true);
                // // println!("meet result fine");
                // // std::mem::forget(left_abs);
                // // meet_ptr
                // left_abs.meet(man, &**right);
                // let meet_ptr = left_abs.elina_abstract0;
                // println!("meet result fine");
                // std::mem::forget(left_abs);
                // meet_ptr
            },
            Binop(Or, left, right) => {
                if destructive {
                    let left_res = left.meet_internal(man, other, false);
                    let right_res = right.meet_internal(man, other, true);
                    let join_res = elina_abstract0_join(man.as_manager_ptr(), true_, right_res, left_res);
                    // right has been mutated, left is still alive but useless => free left
                    // dropping to use the bottom guard
                    std::mem::drop(Abstract { elina_abstract0: left_res });
                    join_res
                } else {
                    let left_res = left.meet_internal(man, other, false);
                    let right_res = right.meet_internal(man, other, false);
                    let join_res = elina_abstract0_join(man.as_manager_ptr(), true_, right_res, left_res);
                    // right is the copied final abs, left is still alive but useless => free left
                    // dropping to use the bottom guard
                    std::mem::drop(Abstract { elina_abstract0: left_res });
                    join_res
                }
                // // TODO: Check this block for memory leaks and interactions (mutating and non-mutating)
                // let left_ptr = left.meet_internal(man, other, false);
                // let left_abs = Abstract { elina_abstract0: left_ptr };
                //
                // let right_ptr = right.meet_internal(man, other, destructive);
                // let mut right_abs = Abstract { elina_abstract0: right_ptr };
                //
                // // then do mutating join on result
                // let join_ptr = left_abs.join_internal(man, &right_abs, true);
                // // forget right_abs, because it's internal pointer is the result of this function
                // std::mem::forget(right_abs);
                // join_ptr
            },
            Unop(Not, inner) => inner.negation().meet_internal(man, other, destructive),
        }
    }
}

impl Meetable for [&Tcons] {
    unsafe fn meet_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        let mut tcons_arr = elina_tcons0_array_make(0);
        let mut ptrs = self
            .iter()
            .map(|tcons| *tcons.elina_tcons0)
            .collect::<Vec<_>>();
        tcons_arr.p = ptrs.as_mut_ptr();
        tcons_arr.size = self.len() as u64;
        std::mem::forget(ptrs);

        let res_abs_ptr = elina_abstract0_meet_tcons_array(
            man.as_manager_ptr(),
            c_bool_from_bool(destructive),
            other,
            &mut tcons_arr,
        );

        std::mem::drop(Vec::from_raw_parts(tcons_arr.p, self.len(), self.len()));

        res_abs_ptr
    }
}

impl<const N: usize> Meetable for [&Tcons; N] {
    unsafe fn meet_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        self[..].meet_internal(man, other, destructive)
    }
}

impl Meetable for Vec<&Tcons> {
    unsafe fn meet_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        self[..].meet_internal(man, other, destructive)
    }
}

impl Meetable for Tcons {
    unsafe fn meet_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        [self].meet_internal(man, other, destructive)
    }
}

impl Meetable for Abstract {
    unsafe fn meet_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        elina_abstract0_meet(
            man.as_manager_ptr(),
            c_bool_from_bool(destructive),
            other,
            self.elina_abstract0,
        )
    }
}

// TODO: think about what we can actually join with. To join with tcons, do we need to do `(top meet tcons) join x`?
/// An element that [`Abstract`] can join with.
///
/// Outside this trait, member functions are only called by
/// [`Abstract::join`] and [`Abstract::join_copy`].
///
/// # Implementation
///
/// Implementors of this trait must either implement [`Joinable::join_internal`] or both
/// [`Joinable::join_with`] and [`Joinable::join_with_copy`].
pub trait Joinable {
    /// Returns a pointer to the internal join result of `self` with `other`.
    ///
    /// If `destructive` is `true`, `other` is internally mutated and its pointer is returned,
    /// if `destructive` is `false`, `other` is not mutated and a pointer to the newly allocated
    /// internal `elina_abstract0_t` is returned.
    ///
    /// # Panics
    ///
    /// This function is only called by [`Joinable::join_with`] and [`Joinable::join_with_copy`].
    /// It panics when neither `join_internal` nor the two `join_with`'s are implemented.
    ///
    /// # Safety
    ///
    /// This is unsafe, because `other` is internally mutated exactly if `destructive` is true.
    /// The caller must ensure that they have mutable permissions for `other` when calling this
    /// function with `destructive` set to true.
    #[allow(unused_variables)]
    unsafe fn join_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        panic!("Joinable::join_internal's definition must be overridden if you don't provide implementations for join_with and join_with_copy");
    }

    // reverse, self.join_with(other) = elina_abstract0_join...(other, self)
    fn join_with<M: Manager>(&self, man: &M, other: &mut Abstract) {
        unsafe {
            // Mutates "other"
            other.elina_abstract0 = self.join_internal(man, other.elina_abstract0, true);
        }
    }
    fn join_with_copy<M: Manager>(&self, man: &M, other: &Abstract) -> Abstract {
        unsafe {
            let new_abs_ptr = self.join_internal(man, other.elina_abstract0, false);
            Abstract {
                elina_abstract0: new_abs_ptr,
            }
        }
    }
}

impl Joinable for Abstract {
    unsafe fn join_internal<M: Manager>(
        &self,
        man: &M,
        other: *mut elina_abstract0_t,
        destructive: bool,
    ) -> *mut elina_abstract0_t {
        elina_abstract0_join(
            man.as_manager_ptr(),
            c_bool_from_bool(destructive),
            other,
            self.elina_abstract0,
        )
    }
}


/// An interval used by [`Abstract::get_bounds`].
///
/// Either bound may be `+/-infinity`, and both bounds are closed.
#[derive(Clone, Copy)]
pub struct Interval(pub Bound, pub Bound);

impl Debug for Interval {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str("[")?;
        Debug::fmt(&self.0, f)?;
        f.write_str(",")?;
        Debug::fmt(&self.1, f)?;
        f.write_str("]")
    }
}

/// A bound used for [`Interval`].
#[derive(Clone, Copy)]
pub enum Bound {
    PosInfinity,
    NegInfinity,
    Num(i64),
}

impl Debug for Bound {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Bound::PosInfinity => f.write_str("+oo"),
            Bound::NegInfinity => f.write_str("-oo"),
            Bound::Num(n) => Debug::fmt(n, f),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::ast::*;

    fn gen_env(num: usize) -> Environment {
        Environment::new((0..num).into_iter().map(|i| format!("x{}", i)).collect())
    }

    #[test]
    fn abstract_meet_tcons() {
        let env = gen_env(2);
    }
}
