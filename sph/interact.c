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
    int i, j, mbins;
    const int n = s->n;
    particle_t* node_buffer[n];
    particle_t* curr;

    float dx, dy, r2, z;
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
    check_bins(s);

    for (i = 0; i < n; ++i) {
        s->particles[i].rho = 4 * s->mass / M_PI / h2;
        get_neighboring_bins(s, &(s->particles[i]), node_buffer, &mbins);

        for (j = 0; j < mbins; ++j) {
            curr = node_buffer[j];
            while (curr) {
                dx = s->particles[i].x[0] - curr->x[0];
                dy = s->particles[i].x[1] - curr->x[1];
                r2 = dx*dx + dy*dy;
                z = h2 - r2;
                if (z > 0 && r2 != 0)
                    s->particles[i].rho += C*z*z*z;

                curr = curr->next;
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
    const float h    = params->h;
    const float rho0 = params->rho0;
    const float k    = params->k;
    const float mu   = params->mu;
    const float g    = params->g;
    const float mass = state->mass;
    const float h2   = h*h;
    int i, j, mbins;

    // Unpack system state
    int n = state->n;

    // Compute density and color
    compute_density(state, params);

    // Start with gravity and surface forces
    for (i = 0; i < n; ++i) {
        state->particles[i].a[0] = 0;
        state->particles[i].a[1] = -g;
    }

    // Constants for interaction term
    float C0 = mass / M_PI / ( (h2)*(h2) );
    float Cp =  15*k;
    float Cv = -40*mu;

    float x, y, dx, dy, r2;
    particle_t* node_buffer[9];
    particle_t* curr;

    for (i = 0; i < n; ++i) {
        const float rhoi = state->particles[i].rho;
        x = state->particles[i].x[0];
        y = state->particles[i].x[1];

        get_neighboring_bins(state, &(state->particles[i]), node_buffer, &mbins);

        for (j = 0; j < mbins; ++j) {
            curr = node_buffer[j];
            while (curr) {
                dx = x - (curr -> x[0]);
                dy = y - (curr -> x[1]);
                r2 = dx*dx + dy*dy;

                if (r2 > 0 && r2 < h2) {
                    const float rhoj = curr -> rho;
                    float q = sqrt(r2)/h;
                    float u = 1 - q;
                    float w0 = C0 * u/rhoi/rhoj;
                    float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
                    float wv = w0 * Cv;
                    float dvx = (state->particles[i].v[0]) - (curr -> v[0]);
                    float dvy = (state->particles[i].v[1]) - (curr -> v[1]);
                    state->particles[i].a[0] += (wp*dx + wv*dvx);
                    state->particles[i].a[1] += (wp*dy + wv*dvy);
                }
                curr = curr -> next;
            }
        }
    }
}

