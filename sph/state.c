#include <stdlib.h>
#include "state.h"
#include "math.h"


sim_state_t* alloc_state(int n)
{
    sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
    s->n     = n;
    s->binw  = (int) floor(sqrt(n));
    s->nbins = s->binw * s->binw;
    s->rho   = (float*) calloc(  n, sizeof(float));
    s->x     = (float*) calloc(2*n, sizeof(float));
    s->vh    = (float*) calloc(2*n, sizeof(float));
    s->v     = (float*) calloc(2*n, sizeof(float));
    s->a     = (float*) calloc(2*n, sizeof(float));
    s->bins  = (node**) calloc(s->nbins, sizeof(node*));

    return s;
}

void clear_bins(sim_state_t* s) {
    int i;
    node *curr, *temp;

    for (i = 0; i < s -> nbins; ++i) {
        curr = s->bins[i];

        while (curr != NULL) {
            temp =  curr->next;
            free(curr);
            curr = temp;
        }
    }
}

void free_state(sim_state_t* s)
{
    free(s->a);
    free(s->v);
    free(s->vh);
    free(s->x);
    free(s->rho);

    clear_bins(s);
    free(s->bins);

    free(s);
}
