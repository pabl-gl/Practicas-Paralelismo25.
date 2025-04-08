#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


int MPI_FlattreeColectiva(const void *sendbuff, void *recbuff, int count, MPI_Datatype type, int root, MPI_Comm comm) {
    int numprocs, rank;
    int received, sum = 0;
    MPI_Status status;
    if (MPI_Comm_size(comm, &numprocs) != MPI_SUCCESS || MPI_Comm_rank(comm, &rank) != MPI_SUCCESS)
        return MPI_ERR_COMM;

    if (rank != root) {
        if (MPI_Send(sendbuff, count, type, root, 1, comm) != MPI_SUCCESS)
            return MPI_ERR_OTHER;
    } else {
        for (int j = 0; j < numprocs; j++) {
            if (j == root) continue;
            if (MPI_Recv(&received, count, type, j, 1, comm, &status) != MPI_SUCCESS)
                return MPI_ERR_OTHER;
            sum += received;
        }
        sum += *(int *) sendbuff; //sum del count del root que no se sumó antes.
        *(int *) recbuff = sum;
    }
    return MPI_SUCCESS;
}

int MPI_BinomialBcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int numprocs, rank, indice;
    if (MPI_Comm_size(comm, &numprocs) != MPI_SUCCESS || MPI_Comm_rank(comm, &rank) != MPI_SUCCESS)
        return MPI_ERR_COMM;
    MPI_Status Status;

    for (int i = 1; pow(2, i - 1) < numprocs; i++) {
        indice = (int) pow(2, i - 1);
        if (rank < indice) {
            int dest = rank + indice;
            if (dest < numprocs) {
                if (MPI_SUCCESS != MPI_Send(buffer, count, datatype, rank + indice, 0, comm))
                    return MPI_ERR_OTHER;
            }
        } else if (rank >= indice && rank < 2 * indice) {
            if (MPI_SUCCESS != MPI_Recv(buffer, count, datatype, rank - indice, 0, comm, &Status))
                return MPI_ERR_OTHER;
        }
    }
    return MPI_SUCCESS;
}


int main(int argc, char *argv[]) {
    int i, done = 0, n, count = 0, rank, numprocs, received;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    while (!done) {
        if (rank == 0) {
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d", &n);
        }
        MPI_BinomialBcast(&n, 1,MPI_INT, 0,MPI_COMM_WORLD);
        if (n == 0) break;
        count = 0;
        for (i = rank + 1; i <= n; i += numprocs) {
            // Get the random numbers between 0 and 1
            x = ((double) rand()) / ((double) RAND_MAX);
            y = ((double) rand()) / ((double) RAND_MAX);

            // Calculate the square root of the squares
            z = sqrt((x * x) + (y * y));

            // Check whether z is within the circle
            if (z <= 1.0)
                count++;
        }
        MPI_FlattreeColectiva(&count, &received, 1,MPI_INT, 0,MPI_COMM_WORLD);
        if (rank == 0) {
            pi = (double) received / (double) n * 4.0;
            printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }
    MPI_Finalize();
    return 0;
}
