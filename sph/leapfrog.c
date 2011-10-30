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

void leapfrog_step(sim_state_t* s, double dt)
{
    int i;
    int n = s->n;

    for (i = 0; i < n; ++i) {
        s->particles[i].vh[0] += s->particles[i].a[0] * dt;
        s->particles[i].vh[1] += s->particles[i].a[1] * dt;
        s->particles[i].v[0]   = s->particles[i].vh[0] + s->particles[i].a[0] * dt / 2;
        s->particles[i].v[0]   = s->particles[i].vh[1] + s->particles[i].a[1] * dt / 2;
        s->particles[i].x[0]  += s->particles[i].vh[0] * dt;
        s->particles[i].x[1]  += s->particles[i].vh[1] * dt;
    }

    /* reflect_bc(s); */
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

    for (i = 0; i < n; ++i) s->particles[i].vh[0] = s->particles[i].v[0] + s->particles[i].a[0] * dt / 2;
    for (i = 0; i < n; ++i) s->particles[i].vh[0] = s->particles[i].v[1] + s->particles[i].a[1] * dt / 2;
    for (i = 0; i < n; ++i) s->particles[i].v[0] += s->particles[i].a[0] * dt;
    for (i = 0; i < n; ++i) s->particles[i].v[1] += s->particles[i].a[1] * dt;
    for (i = 0; i < n; ++i) s->particles[i].x[0] += s->particles[i].vh[0] * dt;
    for (i = 0; i < n; ++i) s->particles[i].x[1] += s->particles[i].vh[1] * dt;

    reflect_bc(s);
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

static void damp_reflect(int which, float barrier, particle_t *offender)
{
    // Coefficient of resitiution
    const float DAMP = 0.75;

    // Ignore degenerate cases
    if (offender->v[which] == 0)
        return;

    // Scale back the distance traveled based on time from collision
    float tbounce = ((offender->x[which])-barrier)/(offender->v[which]);
    offender->x[0] -= (offender->v[0])*(1-DAMP)*tbounce;
    offender->x[1] -= (offender->v[1])*(1-DAMP)*tbounce;

    // Reflect the position and velocity
    offender->x[which]  = 2*barrier-(offender->x[which]);
    offender->v[which]  = -(offender->v[which]);
    offender->vh[which] = -(offender->vh[which]);

    // Damp the velocities
    offender->v[0] *= DAMP;  offender->vh[0] *= DAMP;
    offender->v[1] *= DAMP;  offender->vh[1] *= DAMP;
}


/*@T
 *
 * For each particle, we need to check for reflections on each
 * of the four walls of the computational domain.
 *@c*/
static void reflect_bc(sim_state_t* s)
{
    // Boundaries of the computational domain
    const float XMIN = 0.0;
    const float XMAX = 1.0;
    const float YMIN = 0.0;
    const float YMAX = 1.0;

    int n = s->n;

    for (int i = 0; i < n; ++i) {
        if (s->particles[i].x[0] < XMIN) damp_reflect(0, XMIN, &(s->particles[i]));
        if (s->particles[i].x[0] > XMAX) damp_reflect(0, XMAX, &(s->particles[i]));
        if (s->particles[i].x[1] < YMIN) damp_reflect(1, YMIN, &(s->particles[i]));
        if (s->particles[i].x[1] > YMAX) damp_reflect(1, YMAX, &(s->particles[i]));
    }
}
