#include <stdlib.h>
#include <stdio.h>
#include "state.h"

static void reflect_bc(sim_state_t* s);

/*@T
 * \section{Leapfrog integration}
 *
 * The leapfrog time integration scheme is frequently used in
 * particle simulation algorithms because
 * \begin{itemize}
 * \item It is explicit, which makes it easy to code.
 * \item It is second-order accurate.
 * \item It is {\em symplectic}, which means that it conserves
 *    certain properties of the continuous differential equation
 *    for Hamiltonian systems.  In practice, this means that it
 *    tends to conserve energy where energy is supposed to be
 *    conserved, assuming the time step is short enough for
 *    stability.
 * \end{itemize}
 * Of course, our system is {\em not} Hamiltonian -- viscosity
 * is a form of damping, so the system loses energy.  But we'll
 * stick with the leapfrog integration scheme anyhow.
 *
 * The leapfrog time integration algorithm is named because
 * the velocities are updated on half steps and the positions
 * on integer steps; hence, the two leap over each other.
 * After computing accelerations, one step takes the form
 * \begin{align*}
 *   \bfv^{i+1/2} &= \bfv^{i-1/2} + \bfa^i \Delta t \\
 *   \bfr^{i+1}   &= \bfr^{i}     + \bfv^{i+1/2} \Delta t,
 * \end{align*}
 * This is straightforward enough, except for two minor points.
 * \begin{enumerate}
 * \item
 *   In order to compute the acceleration at time $t$, we need the
 *   velocity at time $t$.  But leapfrog only computes velocities at
 *   half steps!  So we cheat a little: when we compute the half-step
 *   velocity velocity $\bfv^{i+1/2}$ (stored in [[vh]]), we
 *   simultaneously compute an approximate integer step velocity
 *   $\tilde{\bfv}^{i+1}$ (stored in [[v]]) by taking another half
 *   step using the acceleration $\bfa^i$.
 * \item
 *   We don't explicitly represent the boundary by fixed particles,
 *   so we need some way to enforce the boundary conditions.  We take
 *   the simple approach of explicitly reflecting the particles using
 *   the [[reflect_bc]] routine discussed below.
 * \end{enumerate}
 *@c*/

void leapfrog_step(sim_state_t* state, double dt) { /* Seems to be okay */
    int i;
    int n = state -> n;

    for (i = 0; i < n; ++i) {
        particle_node* restrict bin = &(state -> bins[i]);

        bin->vh[0] += bin->a[0] * dt;
        bin->v[0]   = bin->vh[0] + bin->a[0] * dt / 2;
        bin->x[0]  += bin->vh[0] * dt;

        bin->vh[1] += bin->a[1] * dt;
        bin->v[1]   = bin->vh[1] + bin->a[1] * dt / 2;
        bin->x[1]  += bin->vh[1] * dt;
    }

    reflect_bc(state);
}

/*@T
 * At the first step, the leapfrog iteration only has the initial
 * velocities $\bfv^0$, so we need to do something special.
 * \begin{align*}
 *   \bfv^{1/2} &= \bfv^0 + \bfa^0 \Delta t/2 \\
 *   \bfr^{1} &= \bfr^0 + \bfv^{1/2} \Delta t.
 * \end{align*}
 *@c*/


void leapfrog_start(sim_state_t* state, double dt) { /* Looks okay */
    int i;
    int n = state -> n;

    for (i = 0; i < n; ++i) {
        particle_node* restrict bin = &(state -> bins[i]);

        bin->vh[0] = bin->v[0] + bin->a[0] * dt / 2;
        bin->v[0] += bin->a[0] * dt;
        bin->x[0] += bin->vh[0] * dt;

        bin->vh[1] = bin->v[1] + bin->a[1] * dt / 2;
        bin->v[1] += bin->a[1] * dt;
        bin->x[1] += bin->vh[1] * dt;
    }

    reflect_bc(state);
}

/*@T
 *
 * \section{Reflection boundary conditions}
 *
 * Our boundary condition corresponds to hitting an inelastic boundary
 * with a specified coefficient of restitution less than one.  When
 * a particle hits a vertical barrier ([[which = 0]]) or a horizontal
 * barrier ([[which = 1]]), we process it with [[damp_reflect]].
 * This reduces the total distance traveled based on the time since
 * the collision reflected, damps the velocities, and reflects
 * whatever solution components should be reflected.
 *@c*/

static void damp_reflect(int dim, float barrier, particle_node *bin) {
    /* Coefficient of restitution */
    const float DAMP = 0.75;

    /* Ignore degenerate cases */
    if (bin->v[dim] == 0)
      return;

    /* Scale back the distance traveled based on time from collision */
    float tbounce = (bin->x[dim] - barrier)/bin->v[dim];
    bin->x[0] -= bin->v[0] * (1 - DAMP) * tbounce;
    bin->x[1] -= bin->v[1] * (1 - DAMP) * tbounce;

    /* Reflect the position and velocity */
    bin->x[dim]  = 2*barrier - bin->x[dim];
    bin->v[dim]  = -(bin->v[dim]);
    bin->vh[dim] = -(bin->vh[dim]);

    /* Damp the velocities */
    bin->v[0] *= DAMP; bin->vh[0] *= DAMP;
    bin->v[1] *= DAMP; bin->vh[1] *= DAMP;
}
/*@T
 *
 * For each particle, we need to check for reflections on each
 * of the four walls of the computational domain.
 *@c*/
static void reflect_bc(sim_state_t* state) {
    // Boundaries of the computational domain
    const float MIN[2] = {0.0, 0.0};
    const float MAX[2] = {1.0, 1.0};

    int i;
    int n = state -> n;

    for (i = 0; i < n; ++i) {
        particle_node* restrict bin = &(state -> bins[i]);
        if (bin->x[0] < MIN[0]) damp_reflect(0, MIN[0], bin);
        if (bin->x[0] > MAX[0]) damp_reflect(0, MAX[0], bin);
        if (bin->x[1] < MIN[1]) damp_reflect(1, MIN[1], bin);
        if (bin->x[1] > MAX[1]) damp_reflect(1, MAX[1], bin);
    }
}
