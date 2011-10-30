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

void compute_density(sim_state_t* s, sim_param_t* params)
{
    int n = s->n;

    float h  = params->h;
    float h2 = h*h;
    float h8 = ( h2*h2 )*( h2*h2 );
    float C  = 4 * s->mass / M_PI / h8;

    for (int i = 0; i < n; ++i)
        s->particles[i].rho = 0;

    for (int i = 0; i < n; ++i) {
        s->particles[i].rho += 4 * s->mass / M_PI / h2;
        for (int j = i+1; j < n; ++j) {
            float dx = s->particles[i].x[0] - s->particles[j].x[0];
            float dy = s->particles[i].x[1] - s->particles[j].x[1];
            float r2 = dx*dx + dy*dy;
            float z  = h2-r2;
            if (z > 0) {
                float rho_ij = C*z*z*z;
                s->particles[i].rho += rho_ij;
                s->particles[j].rho += rho_ij;
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

    // Unpack system state
    int n = state->n;

    // Compute density and color
    compute_density(state, params);

    // Start with gravity and surface forces
    for (int i = 0; i < n; ++i) {
        state->particles[i].a[0] = 0;
        state->particles[i].a[1] = -g;
    }

    // Constants for interaction term
    float C0 = mass / M_PI / ( (h2)*(h2) );
    float Cp =  15*k;
    float Cv = -40*mu;

    // Now compute interaction forces
    for (int i = 0; i < n; ++i) {
        const float rhoi = state->particles[i].rho;
        for (int j = i+1; j < n; ++j) {
            float dx = state->particles[i].x[0] - state->particles[j].x[0];
            float dy = state->particles[i].x[1] - state->particles[j].x[1];

            float r2 = dx*dx + dy*dy;
            if (r2 < h2) {
                const float rhoj = state->particles[j].rho;
                float q = sqrt(r2)/h;
                float u = 1-q;
                float w0 = C0 * u/rhoi/rhoj;
                float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
                float wv = w0 * Cv;
                float dvx = (state->particles[i].v[0]) - (state->particles[j].v[0]);
                float dvy = (state->particles[i].v[1]) - (state->particles[j].v[1]);
                state->particles[i].a[0] = (wp*dx + wv*dvx);
                state->particles[i].a[1] = (wp*dy + wv*dvy);
                state->particles[j].a[0] = (wp*dx + wv*dvx);
                state->particles[j].a[1] = (wp*dy + wv*dvy);
            }
        }
    }
}

