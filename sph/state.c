#include <stdlib.h>
#include "state.h"

int calc_bins(int n) {
    return 100;
}

sim_state_t* alloc_state(int n) {
    sim_state_t* state = (sim_state_t*) calloc(1, sizeof(sim_state_t));

    state -> n         = n;
    state -> nbins     = calc_bins(n);
    state -> bins      = (particle_node*) calloc(n, sizeof(particle_node));

    return state;
}

void free_state(sim_state_t* state) {
    free(state -> bins);
    free(state);
}
