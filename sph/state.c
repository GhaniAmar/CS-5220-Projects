#include <stdlib.h>
#include "state.h"
#include "math.h"
#include <assert.h>
#include <stdio.h>

#define SQ(x) ((x) * (x))
#define SQDIST(x1,y1,x2,y2) (((x1) - (x2)) * ((x1) - (x2)) + ((y1) - (y2)) * ((y1) - (y2)))

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

void get_neighboring_bins(sim_state_t* state, particle_t* particle, particle_t** node_buffer, int* mbins) {
    int k, width, row, column;
    float x, y, h, h2, binwidth, N, S, E, W, NW, NE, SE, SW;
    k = 0;

    x = particle -> x[0];
    y = particle -> x[1];
    h = state -> h;
    h2 = SQ(h);
    binwidth = 1/((float) state->nbinwidth);

    width = state->nbinwidth;
    column = (int) ((particle->x[0]) * width);
    row = (int) ((particle->x[1]) * width);

    N = y - row * binwidth;
    S = (row + 1) * binwidth - y;
    E = (column + 1) * binwidth - x;
    W = x - column * binwidth;

    node_buffer[k++] = state->bins[column * width + row];

    if (row > 0 && N < h) {
        node_buffer[k++] = state->bins[column * width + (row - 1)];

        NW = SQDIST(x, y, column * binwidth, row * binwidth);
        if (column > 0 && NW < h2)
            node_buffer[k++] = state->bins[(column - 1) * width + (row - 1)];

        NE = SQDIST(x, y, (column + 1) * binwidth, row * binwidth);
        if (column < width - 1 && NE < h2)
            node_buffer[k++] = state->bins[(column + 1) * width + (row - 1)];
    }

    if (row < width - 1 && S < h) {
        node_buffer[k++] = state->bins[column * width + (row + 1)];

        SW = SQDIST(x, y, column * binwidth, (row + 1) * binwidth);
        if (column > 0 && SW < h2)
            node_buffer[k++] = state->bins[(column - 1) * width + (row + 1)];

        SE = SQDIST(x, y, (column + 1) * binwidth, (row + 1) * binwidth);
        if (column < width - 1 && SE < h2)
            node_buffer[k++] = state->bins[(column + 1) * width + (row + 1)];
    }

    if (column > 0 && W < h)
        node_buffer[k++] = state->bins[(column - 1) * width + row];

    if (column < width - 1 && E < h)
        node_buffer[k++] = state->bins[(column + 1) * width + row];

    *mbins = k;
}

/* Sort particle array to improve bin locality */
void bucket_sort(sim_state_t* state) {
    int i, j;
    particle_t* curr;
    particle_t* new_particles = (particle_t*) calloc(state->n, sizeof(particle_t));

    /* Put each particle into its proper bin */
    clear_bins(state);
    for (i = 0; i < state->n; ++i)
        add_to_bin(state, &(state->particles[i]));

    for (i = j = 0; i < state->nbins; ++i) {
        curr = state->bins[i];
        while (curr) {
            new_particles[j++] = *curr;
            curr = curr -> next;
        }
    }

    state->particles = new_particles;
    free(state->particles);
}

void free_state(sim_state_t* s)
{
    free(s->particles);
    free(s->bins);

    free(s);
}
