void relax(double **u, double **rhs, int n)
/*
  Red-black Gauss-Seidel relaxation for model problem. Updates the current value of the solution
  u[1..n][1..n], using the right-hand side function rhs[1..n][1..n].
*/
{
  int i,ipass,isw,j,jsw=1;
  double h,h2;
  h=1.0/(n-1);
  h2=h*h;
  /* Red and black sweeps.*/
  /* jsw and isw toggle between 1 and 2 and
     determine starting row in each column
     for given pass 
  */
  for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) { 
    isw=jsw;
    for (j=2;j<n;j++,isw=3-isw)
      /*Gauss-Seidel formula.*/
      for (i=isw+1;i<n;i+=2) 
	u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]
		      +u[i][j-1]-h2*rhs[i][j]);
  }
}
