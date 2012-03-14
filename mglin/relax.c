#include "nrutil.h"

void relax(double ***u, double ***rhs, int n, double C, int mode, int num_pass)
/*
// USAGE
// Mode == 1 Jacobi
// Mode == 2 Gauss-Seidel
// Mode == 3 Gauss-Seidel Red/Black
*/
{
    double cc;
    cc = 1/(6*C+1);

	if( mode == 1) {
		/*
		Jacobi
		*/
		double ***u_old;
		u_old = d3tensor(1,n,1,n,1,n);
		// make copy of u into u_old
		copy(u_old,u,n);
		
		int i,j,k,ipass;
		for (ipass=1;ipass<=num_pass;ipass++) {
			for (k=2;k<n;k++)
				for (j=2;j<n;j++)
					for (i=2;i<n;i++)  
						u[i][j][k] = C*cc*(u_old[i+1][j][k]+u_old[i-1][j][k]+u_old[i][j+1][k]+u_old[i][j-1][k]+u_old[i][j][k-1]+u_old[i][j][k+1] + cc*rhs[i][j][k]);
		}
		free_d3tensor(u_old,1,n,1,n,1,n);
	}
	else if(mode == 2) {
		/*
		Gauss-Seidel 
		*/
		int i,j,k,ipass;
		for (ipass=1;ipass<=num_pass;ipass++) {
			for (k=2;k<n;k++)
				for (j=2;j<n;j++)
					for (i=2;i<n;i++)  
						u[i][j][k] = C*cc*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1] + cc*rhs[i][j][k]);
		}
	}
	else {
		/*
		Red-black Gauss-Seidel
		*/
		int i,ipass,isw,j,jsw=1,k,ksw=1;
		/* jsw and isw toggle between 1 and 2 and
		determine starting row in each column
		for given pass */
		for (ipass=1;ipass<=num_pass;ipass++,jsw=3-jsw) { 
			isw=jsw;
			isw=3-isw;
			for (k=1;k<n;k++, isw=3-isw)
				for (j=2;j<n;j++,isw=3-isw)
					/*Gauss-Seidel formula.*/
					for (i=isw+1;i<n;i+=2)  
						u[i][j][k] = C*cc*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1] + cc*rhs[i][j][k]);
		}
	}
	return;
}
