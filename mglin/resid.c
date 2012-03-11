void resid(double ***res, double ***u, double ***rhs, int n,double C)
/*
Returns minus the residual for the model problem. Input quantities are u[1..n][1..n] and
rhs[1..n][1..n], while res[1..n][1..n] is returned.
*/
{
  int i,j,k;
  double h,h3i;
  h=1.0/(n-1);
  h3i=1.0/(h*h*h);
  
  /* Interior points.*/
  for (k=2;j<n;k++) 
  	 for (j=2;j<n;j++) 
    	for (i=2;i<n;i++)
      	    res[i][j][k] = -h3i*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+
                    u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1]-(6.0*u[i][j][k]))+rhs[i][j][k];
      	
  /* Boundary points.*/
  for (i=1;i<=n;i++) 
    res[i][1][1]=res[i][n][n]=res[1][i][1]=res[n][i][n]=res[1][1][i]=res[n][n][i]=0.0;
}
