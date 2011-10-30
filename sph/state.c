#include <stdlib.h>
#include "state.h"
#include "math.h"
#include <assert.h>
#include <stdio.h>

sim_state_t* alloc_state(int n, sim_param_t* params)
{
    int i;

    int nbin_upper_bound = floor(1.0 / (2 * params->h)); /* Bins need to be at least 2h wide */
    int nbin_estimate    = floor(sqrt(n));               /* Want at least one particle per bin, on average */

    sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
    s->n           = n;
    s->h           = params->h;
    s->nbinwidth   = nbin_estimate < nbin_upper_bound ? nbin_estimate : nbin_upper_bound;
    s->nbins       = s->nbinwidth * s->nbinwidth;
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

    return s;
}

void clear_bins(sim_state_t* state) {
    int i;

    for (i = 0; i < state->nbins; ++i) {
        state->bins[i] = NULL;
    }
}

void add_to_bin(sim_state_t* state, particle_t* particle) {
    int width = state->nbinwidth;
    int column = (particle->x[0]) * width;
    int row = (particle->x[1]) * width;

    particle->next = state->bins[column * width + row];
    state->bins[column * width + row] = particle;
}

#define SQ(x) ((x) * (x))
/* node_buffer is a 4-array of pointers to bin heads */
/* If less than four are needed, the rest are NULL   */
void get_neighboring_bins(sim_state_t* state, particle_t* particle, particle_t** node_buffer) {
    int i, j, width, row, column;
    float h, h2, h_lower, h_upper, v_lower, v_upper;
    float nw, ne, se, sw;

    h = state->h;
    h2 = h * h;

    for (i = j = 0; i < 4; ++i)
        node_buffer[i] = NULL;

    width = state->nbinwidth;
    row = (particle->x[1]) * width;
    column = (particle->x[0]) * width;

    node_buffer[j++] = state->bins[column * width + row];

    h_lower = (particle->x[1]) - (row / width);
    h_upper = (row + 1) / width - (particle->x[1]);
    v_lower = (particle->x[0]) - (column / width);
    v_upper = (column + 1) / width - (particle->x[0]);

    if (h_lower < h && row > 0)
        node_buffer[j++] = state->bins[column * width + row - 1];
    else if (h_upper < h && row < (width - 1))
        node_buffer[j++] = state->bins[column * width + row + 1];

    if (v_lower < h && column > 0)
        node_buffer[j++] = state->bins[(column - 1) * width + row];
    else if (v_upper < h && column < (width - 1))
        node_buffer[j++] = state->bins[(column + 1) * width + row];

    nw = SQ((column/width) - particle->x[0]) + SQ((row/width) - particle->x[1]);
    ne = SQ((column + 1)/width - particle->x[0]) + SQ((row/width) - particle->x[1]);
    se = SQ((column/width) - particle->x[0]) + SQ((row + 1)/width - particle->x[1]);
    sw = SQ((column + 1)/width - particle->x[0]) + SQ((row + 1)/width - particle->x[1]);

    if (nw < h2 && column > 0 && row > 0)
        node_buffer[j++] = state->bins[(column - 1) * width + (row - 1)];
    else if (se < h && column < (width - 1) && row << (width - 1))
        node_buffer[j++] = state->bins[(column + 1) * width + (row + 1)];

    if (ne < h2 && column < (width - 1) && row > 0)
        node_buffer[j++] = state->bins[(column + 1) * width + (row - 1)];
    else if (sw < h2 && column > 0 && row < (width - 1))
        node_buffer[j++] = state->bins[(column - 1) * width + (row + 1)];

    assert(j <= 4);
}

void check_bins(sim_state_t* state) {
    int i, count;
    int n = state->n;
    particle_t* curr;

    for (i = count = 0; i < state->nbins; ++i) {
        curr = state->bins[i];
        while (curr) {
            count++;
            curr = curr -> next;
        }
    }

    assert(n == count);
}

void free_state(sim_state_t* s)
{
    free(s->particles);
    free(s->bins);

    free(s);
}
