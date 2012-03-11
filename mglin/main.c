#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
int main(int argc, char **argv){
  int i,j;
  FILE *outfile;
  double **f;
  int n=4097;
  int ncycle=2;
  f = dmatrix(1,n,1,n);
  f[2048][2048]=1.0;
  //  for (i=2;i<n;++i)
  //  for (j=2;j<n;++j)
  //    f[i][j] = 2.0;
  mglin(f,n,ncycle);
  outfile = fopen("soln.dat", "w");
  fwrite(&f[1][1],sizeof(double),n*n,outfile);
  fclose(outfile);
}
