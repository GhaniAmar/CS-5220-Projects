#ifndef STATE_H
#define STATE_H

/*@T
 * \section{System state}
 * As suggested, we are using an array of linked lists as our space partitionig
 * data structure.
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

typedef struct particle_node {
    float rho;                  /* Density at particle */
    float x[2];                 /* Position vector */
    float vh[2];                /* Velocity vector (half step) */
    float v[2];                 /* Velocity vector (full step) */
    float a[2];                 /* Acceleration vector */
    struct particle_node *next;
} particle_node;

typedef struct sim_state_t {
    int n;                      /* Number of particles */
    int nbins;                  /* Number of particle bins */
    float mass;                 /* Particle mass */
    particle_node *bins;        /* Array of particle bins */
} sim_state_t;

sim_state_t* alloc_state(int n);
void free_state(sim_state_t* s);

/*@q*/
#endif /* STATE_H */
