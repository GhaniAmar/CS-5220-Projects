#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mt19937p.h"
#include <mpi.h>

int square(int n,               // Number of nodes
           int* restrict l,     // Partial distance at step s
           int* restrict lnew)  // Partial distance at step s+1
{
    int done = 1;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            int lij = lnew[j*n+i];
            for (int k = 0; k < n; ++k) {
                int lik = l[k*n+i];
                int lkj = l[j*n+k];
                if (lik + lkj < lij) {
                    lij = lik+lkj;
                    done = 0;
                }
            }
            lnew[j*n+i] = lij;
        }
    }
    return done;
}

/* Returns integer indicating whether the local columns are done. */
int mpi_square(int n, int n_columns, int* restrict l, int* restrict lnew) {

    /* Idea is to keep a buffer of max number of columns per thread */
    /* Compute local interactions first */
    /* Set local columns to buffer */
    /* Send current buffer to the right (wrap around if at end) */
    /* Compute contributions */
    /* Repeat until we have computed everyone's contribution */
    return 1;
}

static inline void infinitize(int n, int* l)
{
    for (int i = 0; i < n*n; ++i)
        if (l[i] == 0)
            l[i] = n+1;
}

static inline void deinfinitize(int n, int* l)
{
    for (int i = 0; i < n*n; ++i)
        if (l[i] == n+1)
            l[i] = 0;
}

void shortest_paths(int n, int* restrict l, int* iter)
{
    // Generate l_{ij}^0 from adjacency matrix representation
    infinitize(n, l);
    for (int i = 0; i < n*n; i += n+1)
        l[i] = 0;

    // Repeated squaring until nothing changes
    int* restrict lnew = (int*) calloc(n*n, sizeof(int));
    memcpy(lnew, l, n*n * sizeof(int));
    for (int done = 0; !done; ) {
        done = square(n, l, lnew);
        memcpy(l, lnew, n*n * sizeof(int));
    }
    free(lnew);
    deinfinitize(n, l);
}

int mpi_shortest_paths(int n, int* restrict l, int n_columns) {
    int i, iter, all_done, done;
    int* restrict lnew = (int*) calloc(n * n_columns, sizeof(int));
    memcpy(lnew, l, n * n_columns * sizeof(int));

    /* Infinitize the local columns */
    for (i = 0; i < n_columns * n; ++i)
        if (l[i] == 0)
            l[i] = n + 1;

    /* Iterate until MPI_Allreduce tells us we are done */
    for (iter = all_done = 0; !all_done; ++iter) {
        done = mpi_square(n, n_columns, l, lnew);
        MPI_Allreduce(&done, &all_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        memcpy(l, lnew, n * n_columns * sizeof(int));
    }

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
        int graph[n * n];

        /* Generate random graph */
        for (j = 0; j < n; ++j) {
            for (i = 0; i < n; ++i)
                l[j*n + i] = (genrand(&state) < p);
            l[j*n + j] = 0;
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

/*@T
 * \section{The [[main]] event}
 *@c*/

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
