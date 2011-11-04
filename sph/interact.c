#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "params.h"
#include "state.h"
#include "interact.h"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Density computations}
 *
 * The formula for density is
 * \[
 *   \rho_i = \frac{4m}{\pi h^8} \sum_{j \in N_i} (h^2 - r^2)^3.
 * \]
 * We search for neighbors of node $i$ by checking every particle,
 * which is not very efficient.  We do at least take advange of
 * the symmetry of the update ($i$ contributes to $j$ in the same
 * way that $j$ contributes to $i$).
 *@c*/

void compute_density(sim_state_t* s, sim_param_t* params)
{
    int i, pass;
    const int n = s->n;
    const int nbinwidth = s -> nbinwidth;

    const float h  = params->h;
    const float h2 = h*h;
    const float h8 = ( h2*h2 )*( h2*h2 );
    const float C  = 4 * s->mass / M_PI / h8;

    /* We do rebinning once per iteration here */
    clear_bins(s);
    for (i = 0; i < n; ++i) {
        add_to_bin(s, &(s->particles[i]));
        s->particles[i].rho = 4 * s->mass / M_PI / h2;
    }
    clear_flags(s);

    /* Our approach does every third column in parallel   */
    /* Having two columns between each column makes it    */
    /* so each thread does not step on any other thread's */
    /* toes when updating its neighbors. We store a flag  */
    /* with each particle to take advantage of symmetry.  */
    for (pass = 0; pass < 3; ++pass) {

        #pragma omp parallel for shared(s)
        for (i = pass; i < nbinwidth; i += 3) {
            particle_t *node_buffer[4];
            particle_t *curr_bin, *curr;
            int mbins, j, k;
            float dx, dy, r2, z, rho_ij;

            for (j = 0; j < nbinwidth; ++j) {
                curr_bin = s->bins[i * nbinwidth + j];

                /* Go through each particle in the bin */
                while (curr_bin) {
                    get_neighboring_bins(s, curr_bin, node_buffer, &mbins);

                    /* Go through the bin's neighboring bins */
                    for (k = 0; k < mbins; ++k) {
                        curr = node_buffer[k];

                        /* Go through the selected neighbor bin */
                        while(curr) {

                            /* If we've already seen the second particle, then */
                            /* this particle pair has been done already        */
                            if (!curr->flag) {
                                dx = curr_bin->x[0] - curr->x[0];
                                dy = curr_bin->x[1] - curr->x[1];
                                r2 = dx*dx + dy*dy;
                                z = h2 - r2;

                                if (z > 0 && r2 != 0) {
                                    rho_ij = C*z*z*z;
                                    curr_bin->rho += rho_ij;
                                    curr->rho     += rho_ij;
                                }
                            }

                            curr = curr->next;
                        }
                    }

                    curr_bin->flag = 1;
                    curr_bin = curr_bin -> next;
                }
            }
        }
    }
}


/*@T
 * \subsection{Computing forces}
 *
 * The acceleration is computed by the rule
 * \[
 *   \bfa_i = \frac{1}{\rho_i} \sum_{j \in N_i}
 *     \bff_{ij}^{\mathrm{interact}} + \bfg,
 * \]
 * where the pair interaction formula is as previously described.
 * Like [[compute_density]], the [[compute_accel]] routine takes
 * advantage of the symmetry of the interaction forces
 * ($\bff_{ij}^{\mathrm{interact}} = -\bff_{ji}^{\mathrm{interact}}$)
 * but it does a very expensive brute force search for neighbors.
 *@c*/

void compute_accel(sim_state_t* state, sim_param_t* params)
{
    // Unpack basic parameters
    const float h       = params->h;
    const float rho0    = params->rho0;
    const float k       = params->k;
    const float mu      = params->mu;
    const float g       = params->g;
    const float mass    = state->mass;
    const float h2      = h*h;
    const int nbinwidth = state->nbinwidth;

    int i, pass;

    // Unpack system state
    const int n = state->n;

    // Compute density and color
    compute_density(state, params);
    clear_flags(state);

    // Start with gravity and surface forces
    for (i = 0; i < n; ++i) {
        state->particles[i].a[0] = 0;
        state->particles[i].a[1] = -g;
    }

    // Constants for interaction term
    const float C0 = mass / M_PI / ( (h2)*(h2) );
    const float Cp =  15*k;
    const float Cv = -40*mu;

    for (pass = 0; pass < 3; ++pass) {

        /* Used the same trick as compute_density here for parallelization */
        /* We make three passes over the bins, doing every third bin in    */
        /* parallel. This prevents any blocking between the threads, as    */
        /* cannot access the same particle within a pass.                  */
        #pragma omp parallel for shared(state)
        for (i = pass; i < nbinwidth; i += 3) {
            float x, y, dx, dy, r2, rhoi, rhoj, q, u, w0, wp, wv, dvx, dvy, fx, fy;
            particle_t *node_buffer[4];
            particle_t *curr, *curr_bin;
            int j, l, mbins;

            for (j = 0; j < nbinwidth; ++j) {
                curr_bin = state->bins[i * nbinwidth + j];

                while (curr_bin) {
                    x = curr_bin->x[0];
                    y = curr_bin->x[1];
                    rhoi = curr_bin->rho;
                    get_neighboring_bins(state, curr_bin, node_buffer, &mbins);

                    for (l = 0; l < mbins; ++l) {
                        curr = node_buffer[l];

                        while(curr) {
                            if (!curr->flag) {
                                dx = x - curr->x[0];
                                dy = y - curr->x[1];
                                r2 = dx*dx + dy*dy;

                                if (r2 > 0 && r2 < h2) {
                                    rhoj = curr->rho;
                                    q = sqrt(r2)/h;
                                    u = 1 - q;
                                    w0 = C0 * u/rhoi/rhoj;
                                    wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
                                    wv = w0 * Cv;
                                    dvx = curr_bin->v[0] - curr->v[0];
                                    dvy = curr_bin->v[1] - curr->v[1];
                                    fx = wp*dx + wv*dvx;
                                    fy = wp*dy + wv*dvy;
                                    curr_bin->a[0] += fx;
                                    curr_bin->a[1] += fy;
                                    curr->a[0]     -= fx;
                                    curr->a[1]     -= fy;
                                }
                            }

                            curr = curr->next;
                        }
                    }

                    curr_bin->flag = 1;
                    curr_bin = curr_bin->next;
                }
            }
        }
    }
}

