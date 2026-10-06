#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("Outer MPI: Rank %d of %d \n", rank, size);

    // Each rank spawns an inner MPI job
    char command[256];
    snprintf(command, sizeof(command), "mpiexec -n 2 ./mpi_hello");
    int status = system(command);

    if (status != 0) {
        fprintf(stderr, "Inner srun failed for rank %d\n", rank);
    }

    MPI_Finalize();
    return 0;
}