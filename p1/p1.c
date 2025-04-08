#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


int main(int argc, char *argv[]){
    int i, done = 0, n, count, rank,numprocs,received;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;
    MPI_Status Status;


    MPI_Init (&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int root=0;


    while (!done){   
        if (rank==root){
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d",&n);
            for (int j = 1; j < numprocs; j++)
                MPI_Send(&n, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
        }
        else 
            MPI_Recv(&n,1,MPI_INT,root, 0,MPI_COMM_WORLD,&Status);

        count = 0;  
        if (n==0) break;

        for (i = rank+1 ; i <= n; i+=numprocs) {
            // Get the random numbers between 0 and 1
	    x = ((double) rand()) / ((double) RAND_MAX);
	    y = ((double) rand()) / ((double) RAND_MAX);

	    // Calculate the square root of the squares
	    z = sqrt((x*x)+(y*y));

	    // Check whether z is within the circle
	    if(z <= 1.0)
            count++;
        }

        if (rank!=root){
            MPI_Send(&count,1,MPI_INT,root,1,MPI_COMM_WORLD);
        }
        else{
            for (i=1;i<numprocs;i++){
                MPI_Recv(&received,1,MPI_INT,i,1, MPI_COMM_WORLD,&Status);
                count+=received;
            }
            pi = ((double) count/(double) n)*4.0;
            printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }  
    MPI_Finalize();
    return 0;
}


