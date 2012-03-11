void addint(double **uf, double **uc, double **res, int nf)
/*
  Does coarse-to-fine interpolation and adds result to uf. nf is the fine-grid dimension. The
  coarse-grid solution is input as uc[1..nc][1..nc], where nc = nf/2+1. The fine-grid solution
  is returned in uf[1..nf][1..nf]. res[1..nf][1..nf] is used for temporary storage.
*/
{
  void interp(double **uf, double **uc, int nf);
  int i,j;
  interp(res,uc,nf);
  for (j=1;j<=nf;j++)
    for (i=1;i<=nf;i++)
      uf[i][j] += res[i][j];
}
