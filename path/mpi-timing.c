#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mt19937p.h"
#include <mpi.h>

/* We use the ring sharing pattern to pass columns of the original matrix around.        */
/* The columns are passed to the left (n-1) times until everyone has seen everything.    */
/* Returns whether this thread is done. In the iteration loop, we reduce the done flags. */
int mpi_square(int n,
               int* restrict l,
               int* restrict lnew,
               int* restrict in_buffer,
               int* restrict out_buffer,
               double* comp_time,
               double* comm_time) {

    int i, j, k, lij, lik, lkj, rank, size, rank_incr;
    int local_columns, column, column_start, column_end, done;
    int donor, recipient, prev_column, in_columns, out_columns;
    double start = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    start = MPI_Wtime();

    column_start  = (rank * n) / size;
    column_end    = ((rank + 1) * n) / size;
    local_columns = column_end - column_start;

    /* Donor is thead we're receiving from, recipient is thread we're giving to */
    donor = (rank + 1) % size;
    recipient = (size + rank - 1) % size;

    done = 1;

    /* Do local computations (same indexing shenanigans as detailed below) */
    for (j = 0; j < local_columns; ++j) {
        for (i = 0; i < n; ++i) {
            lij = lnew[j*n + i];

            for (k = column_start; k < column_end; ++k) {
                lik = l[(k - column_start)*n + i];
                lkj = l[j*n + k];


                if (lik + lkj < lij) {
                    lij = lik + lkj;
                    done = 0;
                }
            }

            lnew[j*n + i] = lij;
        }
    }

    /* Copy our columns into the in_buffer */
    memcpy(in_buffer, l, n * local_columns * sizeof(int));

    (*comp_time) += MPI_Wtime() - start;

    for (rank_incr = 1; rank_incr < size; ++rank_incr) {

        start = MPI_Wtime();

        /* Compute column currently inside in_buffer */
        prev_column = (rank + rank_incr - 1) % size;
        out_columns = ((prev_column + 1) * n) / size - (prev_column * n) / size;

        /* Move in_buffer to out_buffer */
        memcpy(out_buffer, in_buffer, n * out_columns * sizeof(int));

        /* Compute number of columns expecting from the right */
        column       = (rank + rank_incr) % size;
        column_start = (column * n) / size;
        column_end   = ((column + 1) * n) / size;
        in_columns   = column_end - column_start;

        (*comp_time) += MPI_Wtime() - start;

        start = MPI_Wtime();

        /* Send out_buffer to recipient, receive in_buffer from donor */
        MPI_Sendrecv(out_buffer, n * out_columns, MPI_INT, recipient, 0,
                     in_buffer, n * in_columns, MPI_INT, donor, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        (*comm_time) += MPI_Wtime() - start;

        start = MPI_Wtime();

        /* This is a little tricky. The i and j indexes are local indexes in  */
        /* columns (i.e. column number starts at zero for all threads). The k */
        /* index, on the other hand, is a global index and is converted when  */
        /* we access the communication buffers.                               */
        for (j = 0; j < local_columns; ++j) {
            for (i = 0; i < n; ++i) {
                lij = lnew[j*n + i];

                for (k = column_start; k < column_end; ++k) {
                    lik = in_buffer[(k - column_start)*n + i];
                    lkj = l[j*n + k];

                    if (lik + lkj < lij) {
                        lij = lik + lkj;
                        done = 0;
                    }
                }

                lnew[j*n + i] = lij;
            }
        }

        (*comp_time) += MPI_Wtime() - start;
    }

    return done;
}

int mpi_shortest_paths(int n, int* restrict l, int n_columns, int column_start,
                       double* comp_time, double* comm_time) {
    int i, j, iter, all_done, done;
    int* lnew = calloc(n * n_columns, sizeof(int));
    int* in_buffer = calloc(n * (n_columns + 1), sizeof(int));
    int* out_buffer = calloc(n * (n_columns + 1), sizeof(int));

    /* Infinitize the local columns */
    for (j = 0; j < n_columns; ++j)
        for (i = 0; i < n; ++i)
            if (l[j * n + i] == 0 && i != j + column_start)
                l[j * n + i] = n + 1;

    memcpy(lnew, l, n * n_columns * sizeof(int));

    /* Iterate until MPI_Allreduce tells us we are done */
    for (iter = all_done = 0; !all_done; ++iter) {
      done = mpi_square(n, l, lnew, in_buffer, out_buffer, comp_time, comm_time);

        /* MPI_LAND, a realm of MPI_COMM_WORLD */
        MPI_Allreduce(&done, &all_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        memcpy(l, lnew, n * n_columns * sizeof(int));
    }

    free(out_buffer);
    free(in_buffer);
    free(lnew);

    /* Deinfinitize the local columns */
    for (i = 0; i < n_columns * n; ++i)
        if (l[i] == n + 1)
            l[i] = 0;

    return iter;
}

int* mpi_gen_graph(int n, double p, int column_start, int n_columns) {
    int i, j, rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int* l = calloc(n * n_columns, sizeof(int));
    struct mt19937p state;
    sgenrand(10302011UL, &state);

    /* Advance the PRNG to the correct state */
    for (i = 0; i < n * column_start; ++i)
        genrand(&state);

    /* Generate the local columns */
    for (j = 0; j < n_columns; ++j) {
        for (i = 0; i < n; ++i)
            l[j*n + i] = (genrand(&state) < p);
        l[j*n + j + column_start] = 0;
    }

    return l;
}

/* CAVEAT: The return value is only the correct checksum at processor 0 */
int fletcher16p(int* data, int count, int rank, int nproc) {
    int sum[2] = {0, 0};
    int i;

    if (rank != 0)
        MPI_Recv(sum, 2, MPI_INT, rank - 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (i = 0; i < count; ++i) {
        sum[0] = (sum[0] + data[i]) % 255;
        sum[1] = (sum[1] + sum[0]) % 255;
    }

    if (nproc != 1) {
        MPI_Send(sum, 2, MPI_INT, (rank + 1) % nproc,
                 0, MPI_COMM_WORLD);

        if (rank == 0)
            MPI_Recv(sum, 2, MPI_INT, nproc - 1, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return (sum [1] << 8) | sum[0];
}

void write_matrix(const char* fname, int n, int* a)
{
    FILE* fp = fopen(fname, "w+");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", fname);
        exit(-1);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            fprintf(fp, "%d ", a[j*n+i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

const char* usage =
    "path.x -- Parallel all-pairs shortest path on a random graph\n"
    "Flags:\n"
    "  - n -- number of nodes (200)\n"
    "  - p -- probability of including edges (0.05)\n"
    "  - i -- file name where adjacency matrix should be stored (none)\n"
    "  - o -- file name where output matrix should be stored (none)\n";

int main(int argc, char** argv)
{
    int n    = 200;            // Number of nodes
    double p = 0.05;           // Edge probability
    const char* ifname = NULL; // Adjacency matrix file name
    const char* ofname = NULL; // Distance matrix file name
    int rank, size, iter;
    double t0 = 0.;;
    double comp_time, comm_time;
    comp_time = comm_time = 0.;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Option processing
    extern char* optarg;
    const char* optstring = "hn:d:p:o:i:";
    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch (c) {
        case 'h':
            if (rank == 0)
                fprintf(stderr, "%s", usage);
            return -1;
        case 'n': n = atoi(optarg); break;
        case 'p': p = atof(optarg); break;
        case 'o': ofname = optarg;  break;
        case 'i': ifname = optarg;  break;
        }
    }

    const int column_start = (rank * n) / size;         /* Beginning column index in global matrix */
    const int column_end   = ((rank + 1) * n) / size;   /* End column index in global matrix */
    const int n_columns    = column_end - column_start; /* Number of columns this thread owns */

    int fp;
    int* l = mpi_gen_graph(n, p, column_start, n_columns);

    if (rank == 0) {
        if (ifname)
            printf("This version does not support exporting the generated graph.\n");

        t0 = MPI_Wtime();
    }

    iter = mpi_shortest_paths(n, l, n_columns, column_start, &comp_time, &comm_time);
    fp = fletcher16p(l, n * n_columns, rank, size);

    if (rank == 0)
      printf("Computation: %f\nCommunication: %f\nTotal: %f\n",
             comp_time, comm_time, MPI_Wtime() - t0);

    /* Clean up whatever is left */
    free(l);
    MPI_Finalize();

    return 0;
}        
