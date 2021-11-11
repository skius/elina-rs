#include "opt_pk.h"

// bug from https://github.com/eth-sri/ELINA/issues/87
int main() {
    elina_manager_t *pk = opt_pk_manager_alloc(0);

    // Swap what's commented out, and there won't be a segfault
//    elina_abstract0_t* bot = elina_abstract0_bottom(pk, 1, 0);
    elina_abstract0_t* bot = elina_abstract0_top(pk, 1, 0);

    elina_abstract0_t* top = elina_abstract0_top(pk, 1, 0);

    elina_tcons0_array_t arr = elina_tcons0_array_make(1);
    elina_tcons0_t unsat = elina_tcons0_make_unsat();
    arr.p[0] = unsat;

    elina_abstract0_t* top_meet_unsat = elina_abstract0_meet_tcons_array(pk, true, bot, &arr);
    elina_abstract0_t* top_meet_unsat_join = elina_abstract0_join(pk, true, top_meet_unsat, top);
    return 0;
}