#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int MPI_BinomialBcast (  void * buf , int count , MPI_Datatype datatype , // no terminada
    int root , MPI_Comm comm ){
    
    int numprocs, rank,arbol;
    MPI_Status Status;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); // Obtener el número de procesos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Obtener el ID del proceso actual
    if(rank!=root){
        MPI_Recv(&buf,count,MPI_INT,MPI_ANY_SOURCE,1,comm,&Status);
    }
    for (int i = 1; i < numprocs; i++)
    {
        arbol=2^(i-1);
        if(rank<arbol){
            MPI_Send(&buf,count,MPI_INT,(rank+arbol),1,comm);
        }

    }    
    



    }

int MPI_FlattreeColectiva ( void * buff , void * recvbuff , int count ,
    MPI_Datatype datatype , int root , MPI_Comm comm ) {

    int numprocs, rank,received;
    MPI_Status Status;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); // Obtener el número de procesos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Obtener el ID del proceso actual
    if (rank!=root){
        MPI_Send(&buff,count,MPI_INT,0,1,comm);
    }
    else{
        *(int *)recvbuff = *(int *)buff;
        recvbuff=0;
        for (int i=1;i<numprocs;i++){
            MPI_Recv(&received,count,MPI_INT,i,1, comm,&Status);
           *(int*) recvbuff+=received;
        }
    }
    return 0;


    }

int main(int argc, char *argv[]){
    int i, done = 0, n, count, rank,numprocs,received;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;
    MPI_Status Status;

    MPI_Init (&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);


    while (!done){   
        if (rank==0){
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d",&n);
            MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        }
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
        MPI_FlattreeColectiva(&count,&received,1,MPI_INT,0,MPI_COMM_WORLD);
      if(rank==0){
            pi = ((double) received/(double) n)*4.0;
            printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
      }
    }  
    MPI_Finalize();
    return 0;
}


