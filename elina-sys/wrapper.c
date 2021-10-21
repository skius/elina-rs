#include "wrapper.h"
#include "gmp.h"
#include "mpfr.h"

void foreach_linterm_of_linexpr0(elina_linexpr0_t* e, void (*f) (size_t i, elina_dim_t dim, elina_coeff_t* coeff)) {
    size_t i = 0;
    elina_coeff_t* coeff;
    elina_dim_t dim;

    elina_linexpr0_ForeachLinterm(e,i,dim,coeff){
        f(i, dim, coeff);
    }
}