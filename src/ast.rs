use std::borrow::Borrow;
use std::collections::HashMap;
use std::ffi::{CStr, CString};
use std::fmt::{Debug, Formatter};

use std::ops::Deref;
use std::os::raw::c_char;
use std::ptr::{null_mut, slice_from_raw_parts_mut};

pub use elina_sys::{ConsTyp, TexprBinop, TexprUnop};
use elina_sys::{__gmpf_floor, __gmpq_get_str, __gmpq_out_str, __gmpz_export, __gmpz_fdiv_q, __mpz_struct, bool_from_c_bool, c_bool_from_bool, elina_abstract0_add_dimensions, elina_abstract0_assign_texpr, elina_abstract0_bottom, elina_abstract0_bound_dimension, elina_abstract0_copy, elina_abstract0_dimension, elina_abstract0_free, elina_abstract0_is_bottom, elina_abstract0_is_eq, elina_abstract0_is_top, elina_abstract0_join, elina_abstract0_meet, elina_abstract0_meet_tcons_array, elina_abstract0_sat_tcons, elina_abstract0_t, elina_abstract0_to_lincons_array, elina_abstract0_top, elina_abstract0_widening, elina_constyp_t, elina_constyp_t_ELINA_CONS_DISEQ, elina_constyp_t_ELINA_CONS_EQ, elina_constyp_t_ELINA_CONS_SUPEQ, elina_dim_t, elina_dimchange_t, elina_interval_free, elina_lincons0_array_clear, elina_lincons0_array_print, elina_manager_free, elina_manager_t, elina_scalar_free, elina_scalar_infty, elina_scalar_t, elina_tcons0_array_make, elina_tcons0_t, elina_texpr0_binop, elina_texpr0_copy, elina_texpr0_cst_interval_top, elina_texpr0_cst_scalar_int, elina_texpr0_dim, elina_texpr0_free, elina_texpr0_t, elina_texpr0_unop, elina_texpr_op_t, elina_texpr_rdir_t_ELINA_RDIR_ZERO, elina_texpr_rtype_t_ELINA_RTYPE_INT, false_, free, opt_pk_manager_alloc, true_};
use crate::util::lincons0_to_string;

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

    /// Allocates an [`EnvNames`] that refers to the variable names of `self` for use in FFI.
    ///
    /// See `EnvNames` for its invariants.
    pub fn to_env_names(&self) -> EnvNames {
        let mut hash_vec = self.var_to_dim.iter().collect::<Vec<_>>();
        hash_vec.sort_by(|a, b| a.1.cmp(b.1));

        let names = hash_vec.
            iter()
            .map(|(v, _)| CString::new(v.as_str()).unwrap().into_raw())
            .collect::<Vec<_>>();

        EnvNames {
            names,
        }
    }
}

impl Deref for Environment {
    type Target = HashMap<String, elina_dim_t>;

    fn deref(&self) -> &Self::Target {
        &self.var_to_dim
    }
}

/// A type that owns an array of [`CString`].
///
/// It is intended to be used with various FFI `elina` functions that require an array of strings
/// representing named variables.
pub struct EnvNames {
    names: Vec<*mut c_char>,
}

impl EnvNames {
    /// Returns a pointer to its owned array of strings for use with FFI.
    ///
    /// The caller must ensure that `self` outlives the pointer that this function returns,
    /// because otherwise it will be pointing to garbage.
    ///
    /// **Note:** Additionally, the returned array may not be modified, despite it being a
    /// `*mut _`. This is because currently all FFI functions take `*mut c_char`, even though
    /// in most cases they are not modifying the string.
    pub fn as_mut_ptr(&mut self) -> *mut *mut c_char {
        self.names.as_mut_ptr()
    }
}

impl Drop for EnvNames {
    fn drop(&mut self) {
        unsafe {
            // println!("Dropping EnvNames");
            for ptr in &self.names {
                    std::mem::drop(CString::from_raw(*ptr));
            }
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

    /// Returns a `Texpr` representing all (integer) values.
    pub fn top() -> Texpr {
        unsafe {
            Texpr {
                elina_texpr0: elina_texpr0_cst_interval_top(),
            }
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
///
/// This type can be useful for meeting with an if-condition of a language with nested expressions.
/// Its operations on an Abstract are defined through multiple lower-level operations. As such,
/// expect `Hcons` operations (especially if the `Hcons` is nested deeply) to be slower than
/// [`Tcons`].
#[derive(Clone)]
pub enum Hcons {
    /// Wraps a [`Tcons`].
    Leaf(Tcons),
    /// A binary operation on two constraints.
    Binop(HconsBinop, Box<Hcons>, Box<Hcons>),
    /// A unary operation on a constraint.
    Unop(HconsUnop, Box<Hcons>),
    /// A constrait that is always true, i.e. meeting with Hcons::Top is a no-op.
    Top,
}

impl Hcons {
    /// Returns a constraint representing the conjunction of `self` and `right`.
    pub fn and(self, right: Hcons) -> Hcons {
        Hcons::Binop(
            HconsBinop::And,
            Box::new(self),
            Box::new(right)
        )
    }

    /// Returns a constraint representing the disjunction of `self` and `right`.
    pub fn or(self, right: Hcons) -> Hcons {
        Hcons::Binop(
            HconsBinop::Or,
            Box::new(self),
            Box::new(right)
        )
    }

    /// Returns a constraint representing the negation of `self`.
    ///
    /// This method just constructs a constraint representing the negation, it does not perform
    /// and translations into `or`'s and `and`'s. See [`Hcons::negation`] for that.
    pub fn not(self) -> Hcons {
        Hcons::Unop(HconsUnop::Not, Box::new(self))
    }

    /// Returns `tcons` wrapped in an `Hcons`.
    pub fn leaf(tcons: Tcons) -> Hcons {
        Hcons::Leaf(tcons)
    }


    /// Returns a cloned constraint representing negated `self`.
    ///
    /// This method actually turns negations into applicable `or`'s and `and`'s by using
    /// De Morgan's laws.
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
            // TODO: does this make sense?
            Top => Top,
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

    pub fn into_hcons(self) -> Hcons {
        Hcons::Leaf(self)
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

impl Into<Hcons> for Tcons {
    fn into(self) -> Hcons {
        Hcons::Leaf(self)
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
        self.assign_dim(man, env[var.borrow()], texpr);
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
        self.assign_copy_dim(man, env[var.borrow()], texpr)
    }

    /// Assigns dimension `dim` to `texpr` in `self`.
    ///
    /// This function can be used to model mutable variables.
    ///
    /// See the copying counterpart at [`Abstract::assign_copy_dim`].
    pub fn assign_dim<M: Manager>(
        &mut self,
        man: &M,
        dim: u32,
        texpr: &Texpr,
    ) {
        unsafe {
            elina_abstract0_assign_texpr(
                man.as_manager_ptr(),
                true_,
                self.elina_abstract0,
                dim,
                texpr.elina_texpr0,
                std::ptr::null_mut(),
            );
        }
    }

    /// Returns a new `Abstract` representing `self` after dimension `dim` has been assigned `texpr`.
    ///
    /// This function can be used to model mutable variables.
    ///
    /// See the mutating counterpart at [`Abstract::assign_dim`].
    pub fn assign_copy_dim<M: Manager>(&self, man: &M, dim: u32, texpr: &Texpr) -> Abstract {
        unsafe {
            Abstract {
                elina_abstract0: elina_abstract0_assign_texpr(
                    man.as_manager_ptr(),
                    false_,
                    self.elina_abstract0,
                    dim,
                    texpr.elina_texpr0,
                    std::ptr::null_mut(),
                ),
            }
        }
    }

    /// Returns a new `Abstract` representing `self` widened with `other`.
    ///
    /// Specifically, this function applies the widening operator to
    /// `self` and (`self` JOIN `other`).
    pub fn widen_copy<M: Manager>(&self, man: &M, other: &Abstract) -> Abstract {
        unsafe {
            // Widening first requires join
            let tmp = self.join_copy(man, other);

            // I believe there is no mutating widening?
            let widened = elina_abstract0_widening(man.as_manager_ptr(), self.elina_abstract0, tmp.elina_abstract0);
            Abstract {
                elina_abstract0: widened
            }
        }
    }

    /// Adds `n` dimensions after dimension `dim` to `self`.
    pub fn add_dims<M: Manager>(&mut self, man: &M, dim: u32, n: usize) {
        unsafe {
            let mut dim = dim;

            elina_abstract0_add_dimensions(
                man.as_manager_ptr(),
                true_,
                self.elina_abstract0,
                &mut elina_dimchange_t {
                    dim: &mut dim,
                    intdim: n as u64,
                    realdim: 0
                },
                false_,
            );
        }
    }

    /// Returns a new `Abstract` representing `self` after `n` dimensions have been added after dimension `dim`.
    pub fn add_dims_copy<M: Manager>(&self, man: &M, dim: u32, n: usize) -> Abstract {
        unsafe {
            let mut dim = dim;

            Abstract {
                elina_abstract0: elina_abstract0_add_dimensions(
                    man.as_manager_ptr(),
                    false_,
                    self.elina_abstract0,
                    &mut elina_dimchange_t {
                        dim: &mut dim,
                        intdim: n as u64,
                        realdim: 0
                    },
                    false_,
                )
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

    /// Returns the bounds of the [`Texpr`] in `self`.
    pub fn get_bounds_texpr<M: Manager>(&self, man: &M, texpr: &Texpr) -> Interval {
        unsafe {
            let dimensions = elina_abstract0_dimension(man.as_manager_ptr(), self.elina_abstract0);
            let n = dimensions.intdim as u32;

            let mut state = self.add_dims_copy(man, n-1, 1);
            state.assign_dim(man, n, texpr);

            state.get_bounds_dim(man, n)
        }
    }

    /// Returns the bounds of dimension `dim` in `self`.
    pub fn get_bounds_dim<M: Manager>(&self, man: &M, dim: u32) -> Interval {
        unsafe {
            let interval_ptr =
                elina_abstract0_bound_dimension(man.as_manager_ptr(), self.elina_abstract0, dim);

            let inf = (*interval_ptr).inf;
            let inf_bound = match elina_scalar_infty(inf) {
                b if b < 0 => Bound::NegInfinity,
                b if b > 0 => Bound::PosInfinity,
                _ => {
                    let inf = *(*inf).val.mpq;
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

                    // lower bound, hence we ceil.
                    Bound::Num(((result_enum as f64) / (result_denom as f64)).ceil() as i64)
                }
            };

            let sup = (*interval_ptr).sup;
            let sup_bound = match elina_scalar_infty(sup) {
                b if b < 0 => Bound::NegInfinity,
                b if b > 0 => Bound::PosInfinity,
                _ => {
                    let sup = *(*sup).val.mpq;
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

                    // upper bound, hence we floor.
                    Bound::Num(((result_enum as f64) / (result_denom as f64)).floor() as i64)
                }
            };

            elina_interval_free(interval_ptr);

            Interval(inf_bound, sup_bound)
        }
    }

    /// Returns the bounds of variable `var` in `self`.
    pub fn get_bounds<M, S>(&self, man: &M, env: &Environment, var: S) -> Interval
    where
        M: Manager,
        S: Borrow<str>,
    {
        unsafe {
            let dim = env.var_to_dim[var.borrow()];
            self.get_bounds_dim(man, dim)
        }
    }

    /// Returns a `String` representation of `self`.
    pub fn to_string<M: Manager>(&self, man: &M, env: &Environment) -> String {
        if self.is_bottom(man) {
            // <empty>
            return "<bottom>".to_owned();
        } else if self.is_top(man) {
            // <universal>
            return "<top>".to_owned();
        }
        unsafe {
            let lincons_arr =
                elina_abstract0_to_lincons_array(man.as_manager_ptr(), self.elina_abstract0);

            let len = lincons_arr.size as usize;
            let lincons_slice = slice_from_raw_parts_mut(lincons_arr.p, len);

            let mut result = "{".to_owned();
            let mut first = true;
            for lincons in &mut *lincons_slice {
                if first {
                    result.push_str(" ");
                    first = false;
                } else {
                    result.push_str("; ");
                }

                result.push_str(&lincons0_to_string(man, env, lincons))
            }
            result.push_str(" }");

            result
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

impl PartialEq for Abstract {
    fn eq(&self, other: &Self) -> bool {
        unsafe {
            bool_from_c_bool(
                elina_abstract0_is_eq((*self.elina_abstract0).man, self.elina_abstract0, other.elina_abstract0)
            )
        }
    }
}

impl Eq for Abstract {}

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
                // res will be a new copy if `destructive` is true, or the mutated `other`
                let res = left.meet_internal(man, other, destructive);

                // we are mutably meeting with `right`, since `res` is already either
                // a) a new copy if `destructive` is false
                // b) the mutated `other`, which is fine since `destructive` must be true.
                let res = right.meet_internal(man, res, true);
                res
            },
            Binop(Or, left, right) => {
                let left_res = left.meet_internal(man, other, false);

                // We are allocating an Abstract of `other MEETCOPY left`, because this will only serve
                // as join partner for `right`. Additionally, this correctly frees the temporary
                // abstract0 it made at the end of the function (Drop trait).
                let abs_tmp = Abstract { elina_abstract0: left_res };

                // right_res will be our result abstract0, if `destructive` is true we modify `other`,
                // otherwise we allocate a new one.
                let right_res = right.meet_internal(man, other, destructive);

                // We are joining right_res with abs_tmp mutably, since right_res is already either
                // a) a new copy if `destructive` is false
                // b) the mutated `other`, which is fine since `destructive` must be true.
                let join_res = abs_tmp.join_internal(man, right_res,true);
                join_res
            },
            Unop(Not, inner) => inner.negation().meet_internal(man, other, destructive),
            Top => {
                // Top is a no-op, so we just need to decide whether we are creating a new copy or not
                if destructive {
                    other
                } else {
                    elina_abstract0_copy((*other).man, other)
                }
            }
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
        let _env = gen_env(2);
    }
}
