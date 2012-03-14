#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define NUMBER_SWEEPS 2

typedef struct timeval time_type;

void print_here()
{
    printf("In %s(%d)\n",__FILE__,__LINE__); 
}

time_type start()
{
    print_here();
    time_type t;
    gettimeofday(&t, NULL);
    return t;
}

double stop(time_type *start_time)
{
    print_here();
    double microseconds;
    time_type tv2;
    gettimeofday(&tv2, NULL);
    microseconds = (double)(tv2.tv_sec - start_time->tv_sec) * 1000; // sec to ms
    microseconds += (double)(tv2.tv_usec - start_time->tv_usec) / 1000; // ms to us

    return microseconds;
}

void mglin_driver(double ***f, double C, int n, int ncycle, int nsteps, int mode, int num_steps)
{
    print_here();
    double totalElapsed = 0.0;
    int i;
    int mid = n/2 + 1;
    for (i=1; i<=nsteps; i++)
    {
        //time_type t = start(t);

        mglin(f,n,ncycle, C, mode, num_steps);

        printf("spike value (f[mid][mid][mid]) = %.6f\n", f[mid][mid][mid]);

        //totalElapsed = stop(&t);
    }
    
    //double avgPerTimestep = totalElapsed / nsteps;

    //printf("Problem size = %d ncycles = %d execTime = %.4f (ms/timestep)",
    //        n, ncycle, avgPerTimestep);
}

int main(int argc, char **argv){
  int i,j,k;
  FILE *outfile;
  double ***f;

  int n=129;
  int ncycle=2;
  f = d3tensor(1,n,1,n,1,n);
  
  // setup C
  double dt = 0.0001;
  double alpha = 0.001;
  double lx = 1.0;
  double dx = lx/n;
  double C = alpha*dt/(dx*dx*dx);

  //for (i=1;i<n;++i)
  //    for (j=1;j<n;++j)
  //        for (k=1;k<n;++k)
  //            f[i][j][k] = 0.0;
  
  int c = n/2 + 1;
  f[c][c][c]=1.0;
  int nteps = 100;
  printf("got here\n");

  mglin_driver(f,C,n,ncycle,nteps,1,NUMBER_SWEEPS);

  //mglin(f,n,ncycle);
  outfile = fopen("soln.dat", "w");
  fwrite(&f[1][1][1],sizeof(double),n*n*n,outfile);
  fclose(outfile);
}

