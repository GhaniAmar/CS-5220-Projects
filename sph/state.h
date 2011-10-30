#ifndef STATE_H
#define STATE_H

/*@T
 * \section{System state}
 *
 * The [[sim_state_t]] structure holds the information for the current
 * state of the system and of the integration algorithm.  The array
 * [[x]] has length $2n$, with [[ x[2*i+0] ]] and [[ x[2*i+1] ]] representing
 * the $x$ and $y$ coordinates of the particle positions.  The layout
 * for [[v]], [[vh]], and [[a]] is similar, while [[rho]] only has one
 * entry per particle.
 *
 * The [[alloc_state]] and [[free_state]] functions take care of storage
 * for the local simulation state.
 *@c*/

typedef struct particle_t {
    float rho;                  /* Density                */
    float x[2];                 /* Positions              */
    float vh[2];                /* Velocities (half step) */
    float v[2];                 /* Velocities (full step) */
    float a[2];                 /* Acceleration           */
    struct particle_t *next;    /* Next particle in bin   */
} particle_t;

typedef struct sim_state_t {
    int n;                 /* Number of particles        */
    float mass;            /* Particle mass              */

    int nbins;
    particle_t* particles; /* Array of particles         */
    particle_t** bins;     /* Array of bin head pointers */

} sim_state_t;

sim_state_t* alloc_state(int n);
void free_state(sim_state_t* s);

/*@q*/
#endif /* STATE_H */
