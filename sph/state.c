#include <stdlib.h>
#include "state.h"
#include "math.h"

#include <stdio.h>

sim_state_t* alloc_state(int n)
{
    sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
    s->n           = n;
    s->rho         = (float*) calloc(  n, sizeof(float));
    s->x           = (float*) calloc(2*n, sizeof(float));
    s->vh          = (float*) calloc(2*n, sizeof(float));
    s->v           = (float*) calloc(2*n, sizeof(float));
    s->a           = (float*) calloc(2*n, sizeof(float));

    s->nbins       = (int) floor(sqrt(n));
    s->particles   = (particle_t*) calloc(s->n, sizeof(particle_t));
    s->bins        = (particle_t**) calloc(s->nbins, sizeof(particle_t*));

    printf("Created %d by %d particle bins.\n", s->nbins, s->nbins);

    return s;
}

void free_state(sim_state_t* s)
{
    free(s->a);
    free(s->v);
    free(s->vh);
    free(s->x);
    free(s->rho);

    free(s->particles);
    free(s->bins);

    free(s);
}
