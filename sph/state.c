#include <stdlib.h>
#include "state.h"
#include "math.h"

#include <stdio.h>

sim_state_t* alloc_state(int n)
{
    int i;

    sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
    s->n           = n;
    s->nbins       = (int) floor(sqrt(n));
    s->particles   = (particle_t*) calloc(s->n, sizeof(particle_t));
    s->bins        = (particle_t**) calloc(s->nbins, sizeof(particle_t*));

    for (i = 0; i < n; ++i) {
        s->particles[i].rho   = 0;
        s->particles[i].x[0]  = s->particles[i].x[1] = 0;
        s->particles[i].vh[0] = s->particles[i].vh[1] = 0;
        s->particles[i].v[0]  = s->particles[i].v[1] = 0;
        s->particles[i].a[0]  = s->particles[i].a[1] = 0;
        s->particles[i].next  = NULL;
    }

    printf("Created %d by %d particle bins.\n", s->nbins, s->nbins);

    return s;
}

void free_state(sim_state_t* s)
{
    free(s->particles);
    free(s->bins);

    free(s);
}
