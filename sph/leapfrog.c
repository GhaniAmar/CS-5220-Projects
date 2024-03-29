#include <stdlib.h>
#include <stdio.h>
#include "state.h"
#include <omp.h>

static void damp_reflect(int which, float barrier, sim_state_t* s, int i);

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

void leapfrog_step(sim_state_t* s, double dt)
{
    int i;
    const int n = s->n;
    const float eps  = 1e-4;
    const float XMIN = eps;
    const float XMAX = 1.0 - eps;
    const float YMIN = eps;
    const float YMAX = 1.0 - eps;

    clear_bins(s);

    /* Manually inlined the reflect_bc code here to cut down the number of parallel sections */
    #pragma omp parallel for shared(s) private(i) schedule(static)
    for (i = 0; i < n; ++i) {
        s->particles[i].vh[0] += s->particles[i].a[0] * dt;
        s->particles[i].vh[1] += s->particles[i].a[1] * dt;
        s->particles[i].v[0]   = s->particles[i].vh[0] + s->particles[i].a[0] * dt  / 2;
        s->particles[i].v[1]   = s->particles[i].vh[1] + s->particles[i].a[1] * dt  / 2;
        s->particles[i].x[0]  += s->particles[i].vh[0] * dt;
        s->particles[i].x[1]  += s->particles[i].vh[1] * dt;
        if (s->particles[i].x[0] < XMIN) damp_reflect(0, XMIN, s, i);
        if (s->particles[i].x[0] > XMAX) damp_reflect(0, XMAX, s, i);
        if (s->particles[i].x[1] < YMIN) damp_reflect(1, YMIN, s, i);
        if (s->particles[i].x[1] > YMAX) damp_reflect(1, YMAX, s, i);
    }
}

/*@T
 * At the first step, the leapfrog iteration only has the initial
 * velocities $\bfv^0$, so we need to do something special.
 * \begin{align*}
 *   \bfv^{1/2} &= \bfv^0 + \bfa^0 \Delta t/2 \\
 *   \bfr^{1} &= \bfr^0 + \bfv^{1/2} \Delta t.
 * \end{align*}
 *@c*/

void leapfrog_start(sim_state_t* s, double dt)
{
    int i;
    int n = s->n;

    const float eps  = 1e-4;
    const float XMIN = eps;
    const float XMAX = 1.0 - eps;
    const float YMIN = eps;
    const float YMAX = 1.0 - eps;

    clear_bins(s);

    #pragma omp parallel for shared(s) private(i) schedule(static)
    for (i = 0; i < n; ++i) {
        s->particles[i].vh[0] = s->particles[i].v[0] + s->particles[i].a[0] * dt / 2;
        s->particles[i].vh[1] = s->particles[i].v[1] + s->particles[i].a[1] * dt / 2;
        s->particles[i].v[0] += s->particles[i].a[0] * dt;
        s->particles[i].v[1] += s->particles[i].a[1] * dt;
        s->particles[i].x[0] += s->particles[i].vh[0] * dt;
        s->particles[i].x[1] += s->particles[i].vh[1] * dt;
        if (s->particles[i].x[0] < XMIN) damp_reflect(0, XMIN, s, i);
        if (s->particles[i].x[0] > XMAX) damp_reflect(0, XMAX, s, i);
        if (s->particles[i].x[1] < YMIN) damp_reflect(1, YMIN, s, i);
        if (s->particles[i].x[1] > YMAX) damp_reflect(1, YMAX, s, i);
    }
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

static void damp_reflect(int which, float barrier, sim_state_t* s, int i)
{
    particle_t* restrict p = &(s->particles[i]);

    // Coefficient of restitiution
    const float DAMP = 0.75;

    // Ignore degenerate cases
    if (p->v[which] == 0)
        return;

    // Scale back the distance traveled based on time from collision
    float tbounce = (p->x[which] - barrier)/(p->v[which]);
    p->x[0] -= (p->v[0])*(1-DAMP)*tbounce;
    p->x[1] -= (p->v[1])*(1-DAMP)*tbounce;

    // Reflect the position and velocity
    p->x[which]  = 2*barrier - (p->x[which]);
    p->v[which]  = -(p->v[which]);
    p->vh[which] = -(p->vh[which]);

    // Damp the velocities
    p->v[0]  *= DAMP;
    p->v[1]  *= DAMP;
    p->vh[0] *= DAMP;
    p->vh[1] *= DAMP;
}
