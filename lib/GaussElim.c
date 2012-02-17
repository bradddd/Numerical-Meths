#include "nrutil.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void mprint(double **matrix, int m, char *label);
void upper_triangulate(double **A, double *b, int m, int pivotFlag);
void vprint(double *vector, int m, char *label);
void back_sub(double **A, double *x, double *b, int m);

/*
int old_main(int argc, char **argv){
	for (int flag=1 ; flag >= 0 ; flag--){
		double **A, *b, *x;
		int n = 5;
		A = dmatrix(1,n,1,n);
		b = dvector(1,n);
		x = dvector(1,n);

		A[1][1] =  100.1      ; A[1][2]= 200.2       ; A[1][3] =  400.4        ; A[1][4] = 300.3      ; A[1][5] = 200.2;
		A[2][1] =  100.1e+20  ; A[2][2]= 201.2e+20   ; A[2][3] =  300.3e+20    ; A[2][4] = 200.2e+20  ; A[2][5] = 100.1e+20;
		A[3][1] =  .0000004   ; A[3][2]= .0000002    ; A[3][3] =.0000002       ; A[3][4] =.0000002    ; A[3][5] =.0000004;
		A[4][1] =   4         ; A[4][2]= 3           ; A[4][3] =  2            ; A[4][4] = 1          ; A[4][5] = 2;
		A[5][1] =   3.00000003e+26   ; A[5][2]= 3.00000003e+26        ; A[5][3] =  2.00000002e+26         ; A[5][4] = 4.00000004e+26       ; A[5][5] = 1.00000001e+26 ;
				
		b[1] = 3703.7+100.1; b[2] = 2802.8e+20+100.1e+20+5e+20; b[3] = .0000032+.0000004; b[4] = 28+4; b[5] = 41.00000001e+26+3.00000003e+26 ;
				
		mprint(A,n,"A original");
		vprint(b,n,"b original");

		upper_triangulate(A,b,n,flag);
		back_sub(A,x,b,n);

		mprint(A,n,"A");
		vprint(b,n,"b");
		vprint(x,n,"x");
		
	}
  return 0;
}
*/

void mprint(double **matrix, int m, char *label){
  int i, j;
  printf("%s:\n",label);

  for (i = 1; i <= m; ++i){
    for (j = 1; j <= m; ++j){
      printf("%10.2f ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n------------------------\n");
}

void vprint(double *vector, int m, char *label){
  int i;
  printf("%s:\n",label);
  
  for (i = 1; i <= m; ++i){
    printf("%10.2f ", vector[i]);
  }
  
  printf("\n------------------------\n");
}

void upper_triangulate(double **A, double *b, int m, int pivotFlag){
  int i,j,k;
  
  double scale;
	
	for (j=1;j<m;++j){           /* loop over columns */
		
		//printf("In col %d\n",j);
		
		// set pivot index
		i = j+1;
		
		if (pivotFlag != 0 ){
			//printf("With Pivoting\n");
			double max = fabs(A[j][j]); // use abs to find element w/ max magnitude
			int rowWithMax = j;
			
			// determine row with max
			for ( int d = 0 ; d < (m-i+1) ; d++ ) {
				//printf("d = %d\n",d);
				if ( fabs(A[i+d][j]) > max ) {
					rowWithMax = i+d;
					max= fabs(A[i+d][j]); 
				}
			}
			//printf("Max value of %f in row %d\n",max, rowWithMax);
			
			if ( max == 0) {
				continue; // so you don't do anything
			}
			else{
				int piv = j;
				if (piv!=rowWithMax){
					
					//printf("About to swap %d for %d\n",piv,rowWithMax);
					//mprint(A,m,"A");
					
					//swap A[i] with A[rowWithMax]
					for ( int c=1 ; c <= m ; ++c ){
						double tmp;
						tmp = A[piv][c];
						A[piv][c] = A[rowWithMax][c];
						A[rowWithMax][c] = tmp;
					}
					// swap b elements
					double tmp = b[piv];
					b[piv] = b[rowWithMax];
					b[rowWithMax] = tmp;
					
					//printf("After swap\n");
					//vprint(b,n,"b");
					//mprint(A,m,"A");
				}
			}
			//mprint(A,m,"A");
		}
		
		for (i=j+1;i<=m;++i){      /* loop over rows beneath pivot */
			//printf("in row %d\n",i);
						
			scale = A[i][j]/A[j][j];
			
			/* zero out based on pivot */
			for (k=1;k<=m;++k){
				A[i][k] = A[i][k] - A[j][k]*scale;
			}
			b[i] = b[i] - b[j]*scale; /* same for b */
		}
		
		
	}
}

void back_sub(double **A, double *x, double *b, int m){
  int i,j;
  x[m] = b[m]/A[m][m];

  for (i=m-1;i>=1;--i){
    x[i] = b[i];
    for (j=i+1;j<=m;++j){
      x[i] -= A[i][j]*x[j];
    }
    x[i]/=A[i][i];
  }
}
