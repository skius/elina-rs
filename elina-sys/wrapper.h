//#define HAS_APRON
#include "opt_pk.h"
#include <stdio.h>
//#include "elina_tcons0.h"

//#include "apron_wrapper.h"

// custom function to iterate over a linexpr0's linterms:
void foreach_linterm_of_linexpr0(elina_linexpr0_t* e, void (*f) (size_t i, elina_dim_t dim, elina_coeff_t* coeff));