/**
 * This C program solves the heat conduction equation using the Red/Black
 * Gauss-Seidel method in parallel implementing MPI and OpenMP.
 *
 * Author: Ilya Ilyankou <ilya.ilyankou@trincoll.edu>
 * Date: May 1, 2016
 * Version: 0.1
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <mpi.h>
#include <omp.h>

#define idx(ii, jj, maxJ) ((ii) * (maxJ) + (jj)) // return index of a 2D-array
#define odd(ii, jj) (((ii) + (jj)) & 1)          // return 1 if sum is odd
#define abs(a) ((a) > 0.0 ? a : -a)              // absolute value
#define change(a, b) (abs((a) - (b)) / (abs((a)) + SMALL)) // check convergence

char *inp, *res;    // I/O filenames
FILE *f;            // file pointer
int N,              // number of cores for openmp
    ichunk,         // # of distributed chunks in I direction
    jchunk,         // # of distributed chunks in J direction
    maxIter;        // max # of iterations
double eps;         // epsilon value
double SMALL;       // added to denominator to avoid division by 0


// Output the correct usage. In main(), followed by error and exit.
static int write_usage(char *progname) {
    char const usage[]= "Usage: %s -n ncores -inp filename -res filename\n";
 	char const opts_help[] =
 		"\t-help            print this message\n"
        "\t-n ncores        max number of cores for SMP\n"
        "\t-ichunk n        number of chunks in I direction\n"
        "\t-jchunk n        number of chunks in J direction\n"
 		"\t-inp filename    name of input file\n"
 		"\t-res filename    name of output file\n";

 	fprintf(stdout, usage, progname);
 	fputs(opts_help, stdout);

 	return 0;
 }

// Parse arguments from the terminal. Return 0 if all correct,
// write usage and exit if unknown parameter encountered.
int parse_args(int argc, char **argv) {
    int i;

    inp = malloc(256);
    res = malloc(256);

    for (i = 1; i < argc; i++) {
        if (!strcmp("-help", argv[i])) {
            write_usage(argv[0]);
            exit(EXIT_SUCCESS);
        }

        else if (!strcmp( "-n", argv[i]) && (++i < argc)) {
            N = atoi(argv[i]);
        }

        else if (!strcmp( "-ichunk", argv[i]) && (++i < argc)) {
            ichunk = atoi(argv[i]);
        }

        else if (!strcmp( "-jchunk", argv[i]) && (++i < argc)) {
            jchunk = atoi(argv[i]);
        }

        else if (!strcmp( "-inp", argv[i]) && (++i < argc)) {
            inp = argv[i];
        }

        else if (!strcmp( "-res", argv[i]) && (++i < argc)) {
            res = argv[i];
        }

        else {
            fprintf(stdout, "Unrecognized argument \"%s\"\n", argv[i]);
            write_usage(argv[0]);
            return -1;
        }
    }

    return 0;
}


// Read input from file; assign height and width of the grid to maxI and maxJ.
// Return a pointer to an array of doubles that contains the grid from the file.
double* read_input(char *file, int *maxI, int *maxJ) {
    int i, j, trash;
    double skip;
    double *arr;

    f = fopen(file, "r");

    if (!f) {
        printf("Error: Cannot open input file. Aborting.\n");
        exit(EXIT_FAILURE);
    }

    fscanf(f, "%d", maxJ);
    fscanf(f, "%d", maxI);
    fscanf(f, "%lf", &eps);
    fscanf(f, "%d", &maxIter);

    arr = malloc(sizeof(double) * (*maxI) * (*maxJ));

    for (i = 0; i < (*maxI); i++) {
        for (j = 0; j < (*maxJ); j++) {
            fscanf(f, "%lf", &arr[idx(i, j, (*maxJ))]);
        }
    }

    fclose(f);

    return arr;
}

// Write output to the 'res' file in the exact same format as input.
int write_result(double *arr, int maxI, int maxJ, int localI, int localJ) {
    int i, j, jch, ich;
    int a_size = localI * localJ;        // # of elements (doubles) in arr

    f = fopen(res, "w+");

    fprintf(f, "%d\n", maxJ);
    fprintf(f, "%d\n", maxI);
    fprintf(f, "%lf\n", eps);
    fprintf(f, "%d\n", maxIter);

    for (jch = 0; jch < jchunk; jch++) {
        // Pointer to an upper left element of each chunk j-chunk.
        double *p = &arr[a_size * ichunk * jch];

        // Skip first two rows of each chunk for all jchunk > 0, as they are
        // former boundaries, that is, the same as the last two rows of the
        // previous jchunk (to the north).
        for (i = (jch == 0) ? 0 : 2; i < localI; i++) {
            for (ich = 0; ich < ichunk; ich++) {
                // Similarly, skipping first two elements of each row for all
                // ichunks = 1, 2, 3, ... Since they equal to the last two
                // elements of each row of a chunk to the east
                for (j = (ich == 0) ? 0 : 2; j < localJ; j++) {
                    fprintf(f, "%lf ", p[idx(ich, i * localJ + j, a_size)]);
                }
            }
            fprintf(f, "\n");
        }
    }

    fclose(f);
    return 0;
}

// Perform Gauss-Seidel method to calculate the average of each element in
// array 'arr' with maxI number of rows and maxJ number of columns.
// Avoid borders of the array, that is, the first and the last rows, and
// the first and the last columns. Return 1 if desired presision eps
// is achieved, 0 otherwise.
int gauss(double *arr, int maxI, int maxJ) {
    int i, j,
        C,              // index of the original, "central" element
        N, E, S, W;     // indices of North, East, South and West elements

    double old, new;

    int not_precise = 0; // Assume this iteration is precise, that is,
                         // all changes are < eps.

    // Red
    #pragma omp parallel for \
            private(j, old, new, C, N, E, S, W) \
            reduction(+:not_precise)

    for (i = 1; i < maxI - 1; i++) {
        for (j = 1; j < maxJ - 1; j++) {
            if (odd(i, j)) continue;

            C = idx(i, j, maxJ);
            N = idx(i-1, j, maxJ);
            E = idx(i, j+1, maxJ);
            S = idx(i+1, j, maxJ);
            W = idx(i, j-1, maxJ);

            old = arr[C];
            new = (arr[N] + arr[S] + arr[W] + arr[E]) / 4.0;

            if (change(new, old) > eps) {
                not_precise = 1;
            }

            arr[C] = new;
        }
    }

    // Black
    #pragma omp parallel for \
            private(j, old, new, C, N, E, S, W) \
            reduction(+:not_precise)

    for (i = 1; i < maxI - 1; i++) {
        for (j = 1; j < maxJ - 1; j++) {
            if (!odd(i, j)) continue;

            C = idx(i, j, maxJ);
            N = idx(i-1, j, maxJ);
            E = idx(i, j+1, maxJ);
            S = idx(i+1, j, maxJ);
            W = idx(i, j-1, maxJ);

            old = arr[C];
            new = (arr[N] + arr[S] + arr[W] + arr[E]) / 4.0;

            if (change(new, old) > eps) {
                not_precise = 1;
            }

            arr[C] = new;
        }
    }

    // Return 1 if precision is achieved, 0 otherwise.
    return !not_precise;
}


// Add '.[rank]' to the string inp_r.
// That is, for inp_r = "input.txt" and rank = 7, modify inp_r to "input.txt.7".
int local_filename(char *inp_r, int rank) {
    char *r = malloc(10);
    sprintf(r, "%d", rank);

    strcpy(inp_r, inp);
    strcat(inp_r, ".");
    strcat(inp_r, r);

    return 0;
}

// Return 1 if East neighbor of 'rank' exists, 0 otherwise.
int has_east(int rank) {
    if (rank % ichunk == ichunk - 1) {
        return 0;
    }

    return 1;
}

// Return 1 if West neighbor of 'rank' exists, 0 otherwise.
int has_west(int rank) {
    if (rank % ichunk == 0) {
        return 0;
    }
    return 1;
}

// Return 1 if North neighbor of 'rank' exists, 0 otherwise.
int has_north(int rank) {
    if (rank < ichunk) {
        return 0;
    }
    return 1;
}

// Return 1 if South neighbor of 'rank' exists, 0 otherwise.
int has_south(int rank) {
    if (rank >= ichunk * (jchunk - 1)) {
        return 0;
    }

    return 1;
}


int main(int argc, char **argv) {
    int i, j,
        my_rank,
        nsize,          // size of COMM_WORLD (number of ranks)
        maxI,           // # of rows of the main A grid
        maxJ,           // # of columns of the main A grid
        // The following three must be the same for all grids (chunks)!
        localI,         // # of rows for each chunk.
        localJ,         // # of columns for each chunk.
        a_size,         // # of elements in a small grid
        // And this is used to determine precision.
        stop,           // 1 if precision eps for that rank is achieved,
                        // 0 otherwise.
        should_continue;


    // Initialize values beforehand in case user does not specify them.
    maxIter = 0;
    eps = 0.0;
    N = 1;

    if (parse_args (argc, argv)) {
		return EXIT_FAILURE;
    }

    // Create some small value for SMALL to be added to denominator
    // to avoid division by 0.
    SMALL = eps / 1000.0;

    // Initialize MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nsize);

    omp_set_num_threads(N);

    if (my_rank == 0) {
        // Besides applying Gauss method for its chunk (as all other ranks),
        // rank 0 is in charge of gathering precision data from all other chunks
        // after each iteration. If eps precision is achieved, tells other ranks
        // to stop immediately. At the end, gathers processed grids
        // from every rank into 'a_buff' and writes the resulting file.

        // Generate a filename of the form "input.txt.0"
        char *inp_r = malloc(256);
        local_filename(inp_r, my_rank);

        // Read the grid from "input.txt.0" into 'a'
        double *a = read_input(inp_r, &localI, &localJ);
        a_size = localI * localJ;

        // Calculate height and width of the resulting grid
        maxI = (localI - 2) * jchunk + 2;
        maxJ = (localJ - 2) * ichunk + 2;

        // Initialize 'stop' buffer. For each rank, will contain 0 or 1
        // indicating whether or not eps precision was achieved.
        int *stop_buff = malloc(sizeof(int) * nsize);

        // Initialize South and East boundary buffers only,
        // as 0 has no North or West neighbor.
        double *south_buff = malloc(sizeof(double) * localJ);
        double *east_buff = malloc(sizeof(double) * localI);

        // Initialize the main
        double *a_buff = malloc(sizeof(double) * a_size * nsize);

        double start = omp_get_wtime();     // Start timing!

        for (i = 0; i < maxIter; i++) {
            stop = gauss(a, localI, localJ);
            should_continue = 0;

            MPI_Gather(&stop, 1, MPI_INT, stop_buff, 1, MPI_INT,
                        0, MPI_COMM_WORLD);

            for (j = 0; j < nsize; j++) {
                if (stop_buff[j] == 0) {
                    should_continue = 1;
                    break;
                }
            }

            if (has_south(my_rank)) {
                int neighbor = my_rank + ichunk;
                double *send_buff = &a[idx(localI-2, 0, localJ)];

                // Exchange the South boundary
                MPI_Sendrecv(send_buff, localJ, MPI_DOUBLE, neighbor, 0,
                             south_buff, localJ, MPI_DOUBLE, neighbor, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Update the South boundary from received buffer
                for (j = 0; j < localJ; j++) {
                    a[idx(localI - 1, j, localJ)] = south_buff[j];
                }
            }

            if (has_east(my_rank)) {
                int neighbor = my_rank + 1;
                double *send_buff = malloc(sizeof(double) * localI);

                // Generate East buffer to send as the last element of each row.
                for (j = 0; j < localI; j++) {
                    send_buff[j] = a[idx(j, localJ - 2, localJ)];
                }

                // Exchange the East boundary
                MPI_Sendrecv(send_buff, localI, MPI_DOUBLE, neighbor, 0,
                             east_buff, localI, MPI_DOUBLE, neighbor, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                free(send_buff);

                // Update the East boundary from received buffer
                for (j = 0; j < localI; j++) {
                    a[idx(j, localJ - 1, localJ)] = east_buff[j];
                }
            }

            MPI_Bcast(&should_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (should_continue == 0) {
                break;
            }
        }

        // Gather each rank's processed array into a_buff
        MPI_Gather(a, a_size, MPI_DOUBLE, a_buff, a_size, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        // Done!
        write_result(a_buff, maxI, maxJ, localI, localJ);
        printf("Finished in %lf seconds.\n", omp_get_wtime() - start);

        free(south_buff);
        free(east_buff);
        free(a);
    }

    else {
        // All ranks from 1 to (nsize - 1) do the same thing: read the grid from
        // the file and iterate updating their grid using Gauss-Seidel method.
        // After each iteration, send precision data to rank 0 and wait
        // for instructions whether to proceed to next iteration or stop.

        int should_continue;

        // Open file of form 'input.txt.[rank]' and read grid into array 'a'.
        char *inp_r = malloc(256);
        local_filename(inp_r, my_rank);
        double *a = read_input(inp_r, &localI, &localJ);
        a_size = localI * localJ;

        // Initialie buffers for all four boundaries, in spite of some ranks
        // using two or three of these.
        double *south_buff = malloc(sizeof(double) * localJ);
        double *north_buff = malloc(sizeof(double) * localJ);
        double *west_buff = malloc(sizeof(double) * localI);
        double *east_buff = malloc(sizeof(double) * localI);

        for (i = 0; i < maxIter; i++) {
            // Update the grid and send precision data (1 or 0) to rank 0.
            stop = gauss(a, localI, localJ);
            MPI_Gather(&stop, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (has_north(my_rank)) {
                int neighbor = my_rank - ichunk;

                // Exchange the North boundary.
                MPI_Sendrecv(a + localJ, localJ, MPI_DOUBLE, neighbor, 0,
                             north_buff, localJ, MPI_DOUBLE, neighbor, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Update the North boundary from received buffer.
                for (j = 0; j < localJ; j++) {
                    a[j] = north_buff[j];
                }
            }

            if (has_south(my_rank)) {
                int neighbor = my_rank + ichunk;
                double *send_buff = &a[idx(localI-2, 0, localJ)];

                // Exchange the South boundary.
                MPI_Sendrecv(send_buff, localJ, MPI_DOUBLE,
                             neighbor, 0, south_buff, localJ, MPI_DOUBLE,
                             neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Update the South boundary from received buffer.
                for (j = 0; j < localJ; j++) {
                    a[idx(localI - 1, j, localJ)] = south_buff[j];
                }
            }

            if (has_west(my_rank)) {
                int neighbor = my_rank - 1;
                double *send_buff = malloc(sizeof(double) * localI);

                // Generate send buffer as the first element of each row.
                for (j = 0; j < localI; j++) {
                    send_buff[j] = a[idx(j, 1, localJ)];
                }

                // Exchange the West boundary.
                MPI_Sendrecv(send_buff, localI, MPI_DOUBLE, neighbor, 0,
                             west_buff, localI, MPI_DOUBLE, neighbor, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                free(send_buff);

                // Update the West boundary from received buffer.
                for (j = 0; j < localI; j++) {
                    a[idx(j, 0, localJ)] = west_buff[j];
                }
            }

            if (has_east(my_rank)) {
                int neighbor = my_rank + 1;
                double *send_buff = malloc(sizeof(double) * localI);

                // Generate buffer to send as the last element of each row.
                for (j = 0; j < localI; j++) {
                    send_buff[j] = a[idx(j, localJ - 2, localJ)];
                }

                // Exchange the East boundary.
                MPI_Sendrecv(send_buff, localI, MPI_DOUBLE, neighbor, 0,
                             east_buff, localI, MPI_DOUBLE, neighbor, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                free(send_buff);

                // Update the East boundary from received buffer.
                for (j = 0; j < localI; j++) {
                    a[idx(j, localJ - 1, localJ)] = east_buff[j];
                }
            }

            // Wait for response from rank 0 whether to continue iterating.
            MPI_Bcast(&should_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (should_continue == 0) {
                break;
            }

        }

        // Send local grid to rank 0 and voila!
        MPI_Gather(a, a_size, MPI_DOUBLE, NULL, 0, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        free(north_buff);
        free(south_buff);
        free(east_buff);
        free(west_buff);
        free(a);
    }

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
