#include "mg.h"
#include "nrutil.h"
#include<stdio.h>

void mglin(double ***u, int n, int ncycle, double cTerm,int mode, int num_pass){
/*
  Full Multigrid Algorithm for solution of the steady state heat
  equation with forcing.  On input u[1..n][1..n] contains the
  right-hand side ρ, while on output it returns the solution.  The
  dimension n must be of the form 2j + 1 for some integer j. (j is
  actually the number of grid levels used in the solution, called ng
  below.) ncycle is the number of V-cycles to be used at each level.
*/
  unsigned int j,jcycle,jj,jpost,jpre,nf,ng=0,ngrid,nn;
  double C = cTerm;

  /*** setup multigrid jagged arrays ***/
  double ***iu[NGMAX+1];   /* stores solution at each grid level */
  double ***irhs[NGMAX+1]; /* stores rhs at each grid level */
  double ***ires[NGMAX+1]; /* stores residual at each grid level */
  double ***irho[NGMAX+1]; /* stores rhs during intial solution of FMG */
  
  /*** use bitshift to find the number of grid levels, stored in ng ***/
  nn=n;                   
  while (nn >>= 1) ng++;     
  
  /*** some simple input checks ***/
  if (n != 1+(1L << ng)) nrerror("n-1 must be a power of 2 in mglin.");
  if (ng > NGMAX) nrerror("increase NGMAX in mglin.");
  
  /***restrict solution to next coarsest grid (irho[ng-1])***/
  nn=n/2+1;
  ngrid=ng-1;
  irho[ngrid]=d3tensor(1,nn,1,nn,1,nn); 
  rstrct_3d(irho[ngrid],u,nn);/* coarsens rhs (u at this point) to irho on mesh size nn */
  
  // change u in the main ^^^^^^^^^^^^^^
  
  /***continue setting up coarser grids down to coarsest level***/
  while (nn > 3) { 
    nn=nn/2+1; 
    irho[--ngrid]=d3tensor(1,nn,1,nn,1,nn);
    rstrct_3d(irho[ngrid],irho[ngrid+1],nn); 
  }
  
  /***now setup and solve coarsest level iu[1],irhs[1] ***/
  nn=3;
  iu[1]=d3tensor(1,nn,1,nn,1,nn);
  irhs[1]=d3tensor(1,nn,1,nn,1,nn);
  slvsml_3d(iu[1],irho[1],C);          /* solve the small system directly */
  free_d3tensor(irho[1],1,nn,1,nn,1,nn);
  ngrid=ng;                       /* reset ngrid to original size */

  for (j=2;j<=ngrid;j++) {        /* loop over coarse to fine, starting at level 2 */
      printf("at grid level %d\n",j);
      nn=2*nn-1;  // n/2+1 getting the next lvl of gridpoints
      iu[j]=d3tensor(1,nn,1,nn,1,nn);     /* setup grids for lhs,rhs, and residual */
      irhs[j]=d3tensor(1,nn,1,nn,1,nn);
      ires[j]=d3tensor(1,nn,1,nn,1,nn);
      interp(iu[j],iu[j-1],nn);
      /* irho contains rhs except on fine grid where it is in u */
      copy(irhs[j],(j != ngrid ? irho[j] : u),nn); 
      /* v-cycle at current grid level */
      for (jcycle=1;jcycle<=ncycle;jcycle++) {  
          /* nf is # points on finest grid for current v-sweep */
          nf=nn;                                  
          for (jj=j;jj>=2;jj--) {                 
              for (jpre=1;jpre<=NPRE;jpre++)  /* NPRE g-s sweeps on the finest (relatively) scale */
                  relax(iu[jj],irhs[jj],nf,C,mode,num_pass);
              resid(ires[jj],iu[jj],irhs[jj],nf,C); /* compute res on finest scale, store in ires */
              nf=nf/2+1;                        /* next coarsest scale */
              rstrct_3d(irhs[jj-1],ires[jj],nf);  /* restrict residuals as rhs of next coarsest scale */
              fill0(iu[jj-1],nf);              /* set the initial solution guess to zero */
          } 
          slvsml_3d(iu[1],irhs[1],C);                  /* solve the small problem exactly */
          nf=3;                                   /* fine scale now n=3 */
          for (jj=2;jj<=j;jj++) {                 /* work way back up to current finest grid */
              nf=2*nf-1;                            /* next finest scale */
              addint(iu[jj],iu[jj-1],ires[jj],nf);  /* inter error and add to previous solution guess */
              for (jpost=1;jpost<=NPOST;jpost++)    /* do NPOST g-s sweeps */
                  relax(iu[jj],irhs[jj],nf,C,mode,num_pass);
          }
      }
  }

  copy(u,iu[ngrid],n);              /* copy solution into input array (implicitly returned) */

  /*** clean up memory ***/
  for (nn=n,j=ng;j>=2;j--,nn=nn/2+1) {       
      free_d3tensor(ires[j],1,nn,1,nn,1,nn);      
      free_d3tensor(irhs[j],1,nn,1,nn,1,nn);
      free_d3tensor(iu[j],1,nn,1,nn,1,nn);
      if (j != ng) free_d3tensor(irho[j],1,nn,1,nn,1,nn);
  }
  free_d3tensor(irhs[1],1,3,1,3,1,3);
  free_d3tensor(iu[1],1,3,1,3,1,3);
}


/*** copy ain[n][n] into aout[n][n] ***/
void copy(double ***aout, double ***ain, int n)
{
  int i,j,k;
  for (i=1;i<=n;i++)
    for (j=1;j<=n;j++)
    	for (k=1;k<=n;k++)
     		aout[k][j][i]=ain[k][j][i];
}

/*** fill u[n][n] with zeros ***/
void fill0(double ***u, int n)
{
  int i,j,k;
  for (j=1;j<=n;j++)
    for (i=1;i<=n;i++)
    	for (k=1;k<=n;k++)
	      u[i][j][k]=0.0;
}
