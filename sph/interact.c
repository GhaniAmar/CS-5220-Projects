#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "params.h"
#include "state.h"
#include "interact.h"
#include <omp.h>

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
    int i, j, len, k;
    int n = s->n;
    particle_t* node_buffer[n];
    particle_t* curr;

    float dx, dy, r2, z, rho_ij;
    float h  = params->h;
    float h2 = h*h;
    float h8 = ( h2*h2 )*( h2*h2 );
    float C  = 4 * s->mass / M_PI / h8;
	int tid, size, startindex, endindex;
	
	#pragma omp parallel shared(s,params,n,h,h2,h8,C,size) private(i,j,len,k,node_buffer,curr,dx,dy,r2,z,rho_ij,tid)
	{
	tid = omp_get_thread_num();
	size =  omp_get_num_threads();
	#pragma omp for
    for (i = 0; i < n; ++i) {
        s->particles[i].rho = 4 * s->mass / M_PI / h2;
    }
	startindex= ceil((tid*s->nbins)/size);
	endindex= ceil((tid+1)*s->nbins)/size);
	
    for (i = startindex; i < endindex; ++i) {
		particle_t* llist = s->bins[i];
		while(llist) {
			int mbins;
			get_neighboring_bins(s, &(s->llist), node_buffer, &mbins);

			for (j = 0; j < mbins; ++j) {
				curr = node_buffer[j];
				while (curr) {
					dx = sllist.x[0] - curr->x[0];
					dy = llist.x[1] - curr->x[1];
					r2 = dx*dx + dy*dy;
					z = h2 - r2;
					if (z > 0 && r2 != 0)
						s->llist.rho += C*z*z*z;

					curr = curr->next;
				}
			}
			llist = llist->next;
		}
	}
	}
    /* for (i = 0; i < n; ++i) { */
    /*     s->particles[i].rho += 4 * s->mass / M_PI / h2; */
    /*     for (j = i+1; j < n; ++j) { */
    /*         float dx = s->particles[i].x[0] - s->particles[j].x[0]; */
    /*         float dy = s->particles[i].x[1] - s->particles[j].x[1]; */
    /*         float r2 = dx*dx + dy*dy; */
    /*         float z  = h2-r2; */
    /*         if (z > 0) { */
    /*             float rho_ij = C*z*z*z; */
    /*             s->particles[i].rho += rho_ij; */
    /*             s->particles[j].rho += rho_ij; */
    /*         } */
    /*     } */
    /* } */
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
    int i, j, mbins, tid, size, startindex, endindex;

    // Unpack system state
    int n = state->n;

    // Compute density and color
    compute_density(state, params);
	
    // Start with gravity and surface forces
	#pragma omp parallel for shared(state) private(i)
    for (i = 0; i < n; ++i) {
        state->particles[i].a[0] = 0;
        state->particles[i].a[1] = -g;
    }

    // Constants for interaction term
    float C0 = mass / M_PI / ( (h2)*(h2) );
    float Cp =  15*k;
    float Cv = -40*mu;

    float x, y, dx, dy, r2;
    particle_t* node_buffer[4] = {0, 0, 0, 0};
    particle_t* curr;

    #pragma omp parallel shared(state,params,n,h,rho0,k,mu,g,mass,h2,C0,Cp,Cvsize) private(i,j,mbins,node_buffer,curr,dx,dy,r2,x,y,tid,startindex,endindex)
	{
	tid = omp_get_thread_num();
	size =  omp_get_num_threads();
	startindex= ceil((tid*s->nbins)/size);
	endindex= ceil((tid+1)*s->nbins)/size);
	
	
    for (i = startindex; i < endindex; ++i) {
		particle_t* llist = s->bins[i];
		if (i==endindex-state->nbinwidth) {
				#pragma omp barrier
			}
		while(llist) {
			const float rhoi = state->llist.rho;
			x = state->llist.x[0];
			y = state->llist.x[1];

			get_neighboring_bins(state, &(state->llist), node_buffer, &mbins);

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
						float dvx = (state->llist.v[0]) - (curr -> v[0]);
						float dvy = (state->llist.v[1]) - (curr -> v[1]);
						state->llist.a[0] += (wp*dx + wv*dvx);
						state->llist.a[1] += (wp*dy + wv*dvy);
					}
					curr = curr -> next;
				}
			}
			llist= llist -> next;
		}
    }
	}
}

