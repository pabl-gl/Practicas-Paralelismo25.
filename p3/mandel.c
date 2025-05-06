/*The Mandelbrot set is a fractal that is defined as the set of points c
in the complex plane for which the sequence z_{n+1} = z_n^2 + c
with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi/mpi.h> 

#define DEBUG 1

#define          X_RESN  1024  /* x resolution */
#define          Y_RESN  1024  /* y resolution */

/* Boundaries of the mandelbrot set */
#define           X_MIN  -2.0
#define           X_MAX   2.0
#define           Y_MIN  -2.0
#define           Y_MAX   2.0

/* More iterations -> more detailed image & higher computational cost */
#define   maxIterations  1000

typedef struct complextype
{
  float real, imag;
} Compl;

static inline double get_seconds(struct timeval t_ini, struct timeval t_end)
{
  return (t_end.tv_usec - t_ini.tv_usec) / 1E6 +
         (t_end.tv_sec - t_ini.tv_sec);
}

int main ( int argc, char *argv[])
{

  /* Mandelbrot variables */
  int i, j, k, trozo, rank, nprocs;
  Compl   z, c;
  float   lengthsq, temp;
  int *vres, *res[Y_RESN];
  double tcomm,tcomp,tcommf=0.0,tcompf=0.0;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  trozo = Y_RESN / nprocs;
  int *localRes = malloc(trozo * X_RESN * sizeof(int));


  /* Timestamp variables */
  struct timeval  ti, tf,tci;

  /* Allocate result matrix of Y_RESN x X_RESN */
if (rank == 0) {
    vres = (int *) malloc(Y_RESN * X_RESN * sizeof(int));
    if (!vres) {
        fprintf(stderr, "Error allocating memory\n");
        return 1;
    }

    for (i = 0; i < Y_RESN; i++) {
        res[i] = vres + i * X_RESN;
    }
} else {
    vres = NULL;
}

  

  /* Start measuring time */
  gettimeofday(&ti, NULL);

  /* Calculate and draw points */
  for(i=rank*trozo; i < (rank+1)*trozo; i++)
  {
    for(j=0; j < X_RESN; j++)
    {
      z.real = z.imag = 0.0;
      c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
      c.imag = Y_MAX - i * (Y_MAX - Y_MIN)/Y_RESN;
      k = 0;

      do
      {    /* iterate for pixel color */
        temp = z.real*z.real - z.imag*z.imag + c.real;
        z.imag = 2.0*z.real*z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real*z.real+z.imag*z.imag;
        k++;
      } while (lengthsq < 4.0 && k < maxIterations);

      if (k >= maxIterations) localRes[i % trozo * X_RESN + j]= 0;
      else localRes[i % trozo * X_RESN + j] = k;
    }
  }

  gettimeofday(&tci, NULL);

  MPI_Gather(localRes,trozo*X_RESN,MPI_INT,vres,trozo*X_RESN,MPI_INT,0,MPI_COMM_WORLD);

  /* End measuring time */
  gettimeofday(&tf, NULL);

  tcomm = get_seconds(tci,tf);
  tcomp = get_seconds(ti,tf);
  MPI_Reduce(&tcomm,&tcommf,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&tcomp,&tcompf,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  /* Print result out */
  if(rank==0){
  if( DEBUG ) {
    for(i=0;i<Y_RESN;i++) {
      for(j=0;j<X_RESN;j++)
              printf("%3d ", res[i][j]);
      printf("\n");
    }
    printf("Tiempo computacional: %lf\n",tcomp);
    printf("Tiempo comunicación: %lf\n",tcomm);
    free(vres); 
  }
  }
  free(localRes);
  MPI_Finalize();
  return 0;
}
