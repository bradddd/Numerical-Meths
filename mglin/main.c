#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
int main(int argc, char **argv){
  int i,j,k;
  FILE *outfile;
  double ***f;
  int n=129;
  int ncycle=2;
  f = d3tensor(1,n,1,n,1,n);
  //for (i=1;i<n;++i)
  //    for (j=1;j<n;++j)
  //        for (k=1;k<n;++k)
  //            f[i][j][k] = 0.0;
  
  int c = n/2 + 1;
  f[c][c][c]=1.0;
  mglin(f,n,ncycle);
  outfile = fopen("soln.dat", "w");
  fwrite(&f[1][1][1],sizeof(double),n*n*n,outfile);
  fclose(outfile);
}
