void relax(double ***u, double ***rhs, int n, double C)
/*
  Red-black Gauss-Seidel relaxation for model problem. Updates the current value of the solution
  u[1..n][1..n], using the right-hand side function rhs[1..n][1..n].
*/
{
    double cc;
    cc = 1/(6*C+1);

    int i,ipass,isw,j,jsw,k,ksw=1;
    /* Red and black sweeps.*/
    /* jsw and isw toggle between 1 and 2 and
       determine starting row in each column
       for given pass 
       */
    for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) { 
        isw=jsw;
        isw=3-isw;
        for (k=1;k<n;k++, isw=3-isw)
            for (j=2;j<n;j++,isw=3-isw)
                /*Gauss-Seidel formula.*/
                for (i=isw+1;i<n;i+=2)  
                    u[i][j][k] = C*cc*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1] + cc*rhs[i][j][k]);
    }
}
