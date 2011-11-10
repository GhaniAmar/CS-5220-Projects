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
               int* restrict out_buffer) {

    int i, j, k, lij, lik, lkj, rank, size, rank_incr;
    int local_columns, column, column_start, column_end, n_columns, done;
    int donor, recipient, prev_column, out_columns;
   
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    column_start  = (rank * n) / size;
    column_end    = ((rank + 1) * n) / size;
    local_columns = column_end - column_start;

    /* Donor is thead we're receiving from, recipient is thread we're giving to */
    donor = (rank + 1) % n;
    recipient = (n + rank - 1) % n;

    done = 1;

    /* Do local computations (same k indexing shenanigans as detailed below) */
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
    memcpy(in_buffer, l, n * n_columns * sizeof(int));

    for (rank_incr = 1; rank_incr < n; ++rank_incr) {
        /* Move in_buffer to out_buffer */
        prev_column = (rank + rank_incr - 1) % n;
        out_columns = ((prev_column + 1) * n) / size - (prev_column * n) / size;
        memcpy(out_buffer, in_buffer, n * out_columns * sizeof(int));

        /* Compute number of columns expecting from the right */
        column       = (rank + rank_incr) % n;
        column_start = (column * n) / size;        
        column_end   = ((column + 1) * n) / size;  
        n_columns    = column_end - column_start;

        /* Send out_buffer to recipient, receive in_buffer from donor */
        MPI_Sendrecv(out_buffer, n * out_columns, MPI_INT, recipient, recipient,
                     in_buffer, n * n_columns, MPI_INT, donor, rank,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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
    }

    return done;
}

mpi_shortest_paths(int n, int* restrict l, int n_columns) {
    int i, iter, all_done, done;
    int* restrict lnew = calloc(n * n_columns, sizeof(int));
    int* restrict in_buffer = calloc(n * (n_columns + 1), sizeof(int));
    int* restrict out_buffer = calloc(n * (n_columns + 1), sizeof(int));
    memcpy(lnew, l, n * n_columns * sizeof(int));

    /* Infinitize the local columns */
    for (i = 0; i < n_columns * n; ++i)
        if (l[i] == 0)
            l[i] = n + 1;

    /* Iterate until MPI_Allreduce tells us we are done */
    for (iter = all_done = 0; !all_done; ++iter) {
        done = mpi_square(n, l, lnew, in_buffer, out_buffer);
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

int* gen_graph(int n, double p)
{
    int* l = calloc(n*n, sizeof(int));
    struct mt19937p state;
    sgenrand(10302011UL, &state);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i)
            l[j*n+i] = (genrand(&state) < p);
        l[j*n+j] = 0;
    }
    return l;
}

/* We generate the graph on the first thread and then send the pieces out */
/* Could refine even more to interleave generating and sending! */
int* mpi_gen_graph(const int n, const double p, const int n_columns) {
    int i, j, id, rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int* l = calloc(n * n_columns, sizeof(int));

    if (rank == 0) {
        struct mt19937p state;
        sgenrand(10302011UL, &state);
        int graph[n * n];

        /* Generate random graph */
        for (j = 0; j < n; ++j) {
            for (i = 0; i < n; ++i)
                graph[j*n + i] = (genrand(&state) < p);
            graph[j*n + j] = 0;
        }

        /* Copy my columns into local buffer */
        for (j = 0; j < n_columns; ++j)
            for (i = 0; i < n; ++i)
                l[j*n + i] = graph[j*n + i];

        /* Send other columns out to minions */
        for (id = 1; id < rank; ++id) {
            const int id_column_start = (id * n) / size;
            const int id_column_end   = ((id + 1) * n) / size;
            const int id_n_columns    = id_column_end - id_column_start;

            MPI_Send(graph + n * id_column_start, id_n_columns * n,
                     MPI_INT, id, 0, MPI_COMM_WORLD);
        }
    }
    else
        /* Receive columns from head thread */
        MPI_Recv(l, n_columns * n, MPI_INT, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return l;
}

int* mpi_gen_graph2(int n, double p, int column_start, int n_columns) {
    int i, j, id, rank, size;
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
      l[j*n + j] = 0;
    }

    /* Make sure everyone has generated before moving on */
    MPI_Barrier(MPI_COMM_WORLD);

    return l;
}

int fletcher16(int* data, int count)
{
    int sum1 = 0;
    int sum2 = 0;
    for(int index = 0; index < count; ++index) {
          sum1 = (sum1 + data[index]) % 255;
          sum2 = (sum2 + sum1) % 255;
    }
    return (sum2 << 8) | sum1;
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
    double t0 = 0;;

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

    int* l = mpi_gen_graph(n, p, n_columns);

    if (rank == 0) {
        if (ifname)
            printf("This version does not support exporting graphs.\n");

        t0 = MPI_Wtime();
    }

    iter = mpi_shortest_paths(n, l, n_columns);

    /* Need to implement checksum that gathers */

    if (rank == 0) {
        printf("== MPI with %d threads\n", size);
        printf("n:         %d\n", n);
        printf("p:         %g\n", p);
        printf("squares:   %d\n", iter);
        printf("Time:      %g\n", MPI_Wtime() - t0);
        /* printf("Check:     %X\n", fletcher16(l, n*n)); */

        /* Generate output file */
        if (ofname)
            printf("This version does not support exporting graphs.\n");
    }

    /* Clean up whatever is left */
    free(l);
    MPI_Finalize();

    return 0;
}
