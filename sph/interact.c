#include <string.h>
#include <math.h>
#include <assert.h>

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

void compute_density(sim_state_t *state, sim_param_t *params) { /* Looks okay */
    int i, j;
    int n = state -> n;

    float dx, dy, r2, z, rho_ij;
    float h  = params -> h;
    float h2 = h * h;
    float h8 = (h2 * h2) * (h2 * h2);
    float C  = 4 * state->mass / M_PI / h8;

    for (i = 0; i < n; ++i)
        state->bins[i].rho = 0;

    for (i = 0; i < n; ++i) {
        particle_node* restrict bin = &(state->bins[i]);
        bin->rho += 4 * state->mass / M_PI / h2;

        for (j = i + 1; j < n; ++j) {
            particle_node* restrict bin2 = &(state -> bins[j]);
            dx = (bin->x[0]) - (bin2->x[0]);
            dy = (bin->x[1]) - (bin2->x[1]);
            r2 = dx * dx + dy * dy;
            z = h2 - r2;
            if (z > 0) {
                rho_ij = C*z*z*z;
                bin->rho  += rho_ij;
                bin2->rho += rho_ij;
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

void compute_accel(sim_state_t* state, sim_param_t* params) { /* Looks fine */
    int i, j;
    float dx, dy, r2, q, u, w0, wp, wv, dvx, dvy;

    /* Unpack basic parameters */
    const float h    = params->h;
    const float rho0 = params->rho0;
    const float k    = params->k;
    const float mu   = params->mu;
    const float g    = params->g;
    const float mass = state->mass;
    const float h2   = h*h;

    /* Unpack system state */
    particle_node* restrict bins = state -> bins;
    int n = state -> n;

    /* Compute desnity and color */
    compute_density(state, params);

    /* Start with gravity and surface forces */
    for (i = 0; i < n; ++i) {
        bins[i].a[0] = 0;
        bins[i].a[1] = -g;
    }

    /* Constants for interaction term */
    float C0 = mass / M_PI / (h2*h2);
    float Cp =  15*k;
    float Cv = -40*mu;

    /* Now compute interaction forces */
    for (i = 0; i < n; ++i) {
        particle_node* restrict bin = &bins[i];
        const float rhoi = bin->rho;

        for (j = i + 1; j < n; ++j) {
            particle_node* restrict bin2 = &bins[j];
            const float rhoj = bin2->rho;

            dx = (bin->x[0]) - (bin2->x[0]);
            dy = (bin->x[1]) - (bin2->x[1]);
            r2 = dx * dx + dy * dy;

            if (r2 < h2) {
                q = sqrt(r2)/h;
                u = 1 - q;
                w0 = C0 * u / rhoi / rhoj;
                wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
                wv = w0 * Cv;
                dvx = bin->v[0] - bin2->v[0];
                dvy = bin->v[1] - bin2->v[1];
                bin->a[0]  += (wp * dx + wv * dvx);
                bin->a[1]  += (wp * dy + wv * dvy);
                bin2->a[0] -= (wp * dx + wv * dvx);
                bin2->a[1] -= (wp * dy + wv * dvy);
            }
        }
    }
}

