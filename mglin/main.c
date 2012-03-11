#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
int main(int argc, char **argv){
  int i,j;
  FILE *outfile;
  double ***f;
  int n=17;
  int ncycle=2;
  f = d3tensor(1,n,1,n,1,n);
  f[8][8][8]=1.0;
  //  for (i=2;i<n;++i)
  //  for (j=2;j<n;++j)
  //    f[i][j] = 2.0;
  mglin(f,n,ncycle);
  outfile = fopen("soln.dat", "w");
  fwrite(&f[1][1][1],sizeof(double),n*n*n,outfile);
  fclose(outfile);
}
