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

int len(particle_t* particle) {
    int n = 0;
    while (particle) {particle = particle->next; n++;}
    return n;
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

#define SQ(x) ((x) * (x))

/* Returns every particle as its own bin */
/* void get_neighboring_bins(sim_state_t* state, particle_t* particle, particle_t** node_buffer, int* mbins) { */
/*     int i; */
/*     int n = state->n; */

/*     for (i = 0; i < n; ++i) { */
/*         state->particles[i].next = NULL; */
/*         node_buffer[i] = &(state->particles[i]); */
/*     } */

/*     *mbins = state->n; */
/* } */

/* Returns the ~9 bins around the particle */
void get_neighboring_bins(sim_state_t* state, particle_t* particle, particle_t** node_buffer, int* mbins) {
    int k, width, row, column;
    k = 0;

    width = state->nbinwidth;
    column = floor((particle->x[0]) * width);
    row = floor((particle->x[1]) * width);

    node_buffer[k++] = state->bins[column * width + row];

    if (row > 0) {
        node_buffer[k++] = state->bins[column * width + (row - 1)];

        if (column > 0)
            node_buffer[k++] = state->bins[(column - 1) * width + (row - 1)];

        if (column < width - 1)
            node_buffer[k++] = state->bins[(column + 1) * width + (row - 1)];
    }

    if (row < width - 1) {
        node_buffer[k++] = state->bins[column * width + (row + 1)];

        if (column > 0)
            node_buffer[k++] = state->bins[(column - 1) * width + (row + 1)];

        if (column < width - 1)
            node_buffer[k++] = state->bins[(column + 1) * width + (row + 1)];
    }

    if (column > 0)
        node_buffer[k++] = state->bins[(column - 1) * width + row];

    if (column < width - 1)
        node_buffer[k++] = state->bins[(column + 1) * width + row];

    *mbins = k;
}

/* Get bins around the particle (pointers written to node_buffer) */
/* void get_neighboring_bins(sim_state_t* state, particle_t* particle, particle_t** node_buffer, int* mbins) { */
/*     int i, j, width, row, column; */
/*     float h, h2, h_lower, h_upper, v_lower, v_upper; */
/*     float nw, ne, se, sw; */

/*     h = state->h; */
/*     h2 = h * h; */

/*     for (i = j = 0; i < 4; ++i) */
/*         node_buffer[i] = NULL; */

/*     width = state->nbinwidth; */
/*     row = floor((particle->x[1]) * width); */
/*     column = floor((particle->x[0]) * width); */

/*     /\* Add the bin containing the particle *\/ */
/*     node_buffer[j++] = state->bins[column * width + row]; */

/*     /\* Compute distances from the particle to the edges of containing bin *\/ */
/*     v_lower = (particle->x[1]) - (row / width); */
/*     v_upper = (row + 1) / width - (particle->x[1]); */
/*     h_lower = (particle->x[0]) - (column / width); */
/*     h_upper = (column + 1) / width - (particle->x[0]); */

/*     /\* If vertical distances are less than h, add the bin above or below *\/ */
/*     if (v_lower < h && row > 0) */
/*         node_buffer[j++] = state->bins[column * width + row - 1]; */
/*     else if (v_upper < h && row < (width - 1)) */
/*         node_buffer[j++] = state->bins[column * width + row + 1]; */

/*     /\* If the horizontal distances are less than h, add the bin to the left or right *\/ */
/*     if (h_lower < h && column > 0) */
/*         node_buffer[j++] = state->bins[(column - 1) * width + row]; */
/*     else if (h_upper < h && column < (width - 1)) */
/*         node_buffer[j++] = state->bins[(column + 1) * width + row]; */

/*     /\* Compute the distances from the particle to the corners of the bin *\/ */
/*     nw = SQ((column/width) - particle->x[0]) + SQ((row/width) - particle->x[1]); */
/*     ne = SQ((column + 1)/width - particle->x[0]) + SQ((row/width) - particle->x[1]); */
/*     se = SQ((column/width) - particle->x[0]) + SQ((row + 1)/width - particle->x[1]); */
/*     sw = SQ((column + 1)/width - particle->x[0]) + SQ((row + 1)/width - particle->x[1]); */

/*     /\* Add the NW and SE diagonal bins if the distances are less than h *\/ */
/*     if (nw < h2 && column > 0 && row > 0) */
/*         node_buffer[j++] = state->bins[(column - 1) * width + (row - 1)]; */
/*     else if (se < h && column < (width - 1) && row < (width - 1)) */
/*         node_buffer[j++] = state->bins[(column + 1) * width + (row + 1)]; */

/*     /\* Do the same for the NE and NW bins *\/ */
/*     if (ne < h2 && column < (width - 1) && row > 0) */
/*         node_buffer[j++] = state->bins[(column + 1) * width + (row - 1)]; */
/*     else if (sw < h2 && column > 0 && row < (width - 1)) */
/*         node_buffer[j++] = state->bins[(column - 1) * width + (row + 1)]; */

/*     /\* Number of bins added to the buffer is j *\/ */
/*     *mbins = j; */

/*     /\* Sanity check (not much left) *\/ */
/*     assert(j <= 4); */
/* } */

void free_state(sim_state_t* s)
{
    free(s->particles);
    free(s->bins);

    free(s);
}
