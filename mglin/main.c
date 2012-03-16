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

enum RelaxtionMode
{
    RELAX_JACOBI = 1,
    RELAX_GS = 2,
    RELAX_GAUSS_GS_RB = 3,
    RELAX_MAX_VALUE
};

enum DiscMethod
{
    DISC_BE = 1,
    DISC_CN = 2,
    DISC_MAX_VALUE
};

enum InterpType
{
    INTERP_TRILIN = 1,
    INTERP_CUBIC = 2,
    INTERP_VALUE
};

enum RestrictType
{
    RESTRICT_WEIGHTED = 1,
    RESTRICT_DIRECT = 2,
    RESTRICT_MAX_VALUE
};


const char *RelaxtionModeToString(int mode)
{
    if (mode == 1)
        return "JACOBI_RELAXATION";
    else if (mode == 2)
        return "GS_RELAXATION";
    else if (mode == 3)
        return "GS_RB_RELAXATION";
    return "UNKNOWN_RELAX_METHOD";
}

const char *DiscMethodToString(int mode)
{
    if (mode == 1)
        return "DISC_BE";
    else if (mode == 2)
        return "DISC_CN";
    return "UNKNOWN_DISC_METHOD";
}

const char *InterpTypeToString(int mode)
{
    if (mode == 1)
        return "INTERP_TRI";
    else if (mode == 2)
        return "INTERP_CUBIC";
    return "UNKNOWN_INTERP_TYPE";
}

const char *RestrictTypeToString(int mode)
{
    if (mode == 1)
        return "RESTRICT_WEIGHTED";
    else if (mode == 2)
        return "RESTRICT_DIRECT";
    return "UNKNOWN_RESTRICT_TYPE";
}



void mglin_driver(FILE *file, double ***f, double C, int n, int ncycle, int nsteps, int mode, int num_steps, int restrict_mode)
{
    double totalElapsed = 0.0;
    int i;
    int mid = n/2 + 1;

    for (i=1; i<=nsteps; i++)
    {
        time_type t1;
        gettimeofday(&t1, NULL);

        mglin(f,n,ncycle, C, mode, num_steps, restrict_mode);
        //printf("spike value (f[mid][mid][mid]) = %.6f\n", f[mid][mid][mid]);

        double microseconds;
        time_type t2;
        gettimeofday(&t2, NULL);
        microseconds = (double)(t2.tv_sec - t1.tv_sec) * 1000; // sec to ms
        microseconds += (double)(t2.tv_usec - t1.tv_usec) / 1000; // ms to us
        totalElapsed += microseconds;
    }

    double avgPerTimestep = totalElapsed / (double)nsteps;

    printf("Relax(%s) Problem size(%d) ncycles(%d) execTime=%.4f (ms/timestep)\n",
            RelaxtionModeToString(mode), n, ncycle, avgPerTimestep);

    char buffer[256] = {0};
    sprintf(buffer, "%s, %d, %d, %.4f\n",RelaxtionModeToString(mode),n,ncycle,avgPerTimestep );
    fwrite(buffer, 256, 1, file);
}

int main(int argc, char **argv){
  FILE *outfile;
  FILE *timefile;
  double ***f;

  int matrix_size[4] = {33, 65, 129, 257};
  int ncycle=2;
  
  timefile = fopen("FMA_benchmark.csv", "w");
  char buf[128]={0};
  sprintf(buf, "RelaxationType, ProblemSize, Nsteps, Time(us)\n");
  fwrite(buf, 128, 1, timefile);

  int k;
  for (k = 0; k < sizeof(matrix_size)/sizeof(int); k++)
  {

      int n = matrix_size[k];
      // setup C
      double dt = 0.0005;
      double alpha = 0.01;
      double lx = 1.0;
      double dx = lx/n;
      double C = alpha*dt/(dx*dx*dx);

      int i;
      for (i = 1; i < RELAX_MAX_VALUE; i++)
      {
          f = d3tensor(1,n,1,n,1,n);
          int c = n/2 + 1;
          f[c][c][c]=1.0;
          int nsteps = 100;

          printf("Running FMA (n = %d, nsteps=%d, C=%.4f\n", n, nsteps, C);
          mglin_driver(timefile, f,C,n,ncycle,nsteps,i,NUMBER_SWEEPS,1);

          outfile = fopen("soln.dat", "w");
          fwrite(&f[1][1][1],sizeof(double),n*n*n,outfile);
          fclose(outfile);

          free_d3tensor(f,1,n,1,n,1,n);
      }
  }
  fclose(timefile);


}

