//scrap code

// class Matrix{
	
	// void mprint(double **matrix, int m, char *label){
	  // int i, j;
	  // printf("%s:\n",label);

	  // for (i = 1; i <= m; ++i){
		// for (j = 1; j <= m; ++j){
		  // printf("%10.2f ", matrix[i][j]);
		// }
		// printf("\n");
	  // }
	  // printf("\n------------------------\n");
	// }

	// void vprint(double *vector, int m, char *label){
	  // int i;
	  // printf("%s:\n",label);
	  
	  // for (i = 1; i <= m; ++i){
		// printf("%10.2f ", vector[i]);
	  // }
	  
	  // printf("\n------------------------\n");
	// }

	// void upper_triangulate(double **A, double *b, int m, int pivotFlag){
		// int i,j,k;
	  
		// double scale;
		
		// for (j=1;j<m;++j){           /* loop over columns */
			
			// //printf("In col %d\n",j);
			
			// // set pivot index
			// i = j+1;
			
			// if (pivotFlag != 0 ){
				// printf("With Pivoting\n");
				// double max = fabs(A[j][j]); // use abs to find element w/ max magnitude
				// int rowWithMax = j;
				
				// // determine row with max
				// for ( int d = 0 ; d < (m-i+1) ; d++ ) {
					// //printf("d = %d\n",d);
					// if ( fabs(A[i+d][j]) > max ) {
						// rowWithMax = i+d;
						// max= fabs(A[i+d][j]); 
					// }
				// }
				// //printf("Max value of %f in row %d\n",max, rowWithMax);
				
				// if ( max == 0) {
					// continue; // so you don't do anything
				// }
				// else{
					// int piv = j;
					// if (piv!=rowWithMax){
						
						// //printf("About to swap %d for %d\n",piv,rowWithMax);
						// //mprint(A,m,"A");
						
						// //swap A[i] with A[rowWithMax]
						// for ( int c=1 ; c <= m ; ++c ){
							// double tmp;
							// tmp = A[piv][c];
							// A[piv][c] = A[rowWithMax][c];
							// A[rowWithMax][c] = tmp;
						// }
						// // swap b elements
						// double tmp = b[piv];
						// b[piv] = b[rowWithMax];
						// b[rowWithMax] = tmp;
						
						// //printf("After swap\n");
						// //vprint(b,n,"b");
						// //mprint(A,m,"A");
					// }
				// }
				// //mprint(A,m,"A");
			// }
			
			// for (i=j+1;i<=m;++i){      /* loop over rows beneath pivot */
				// //printf("in row %d\n",i);
							
				// scale = A[i][j]/A[j][j];
				
				// /* zero out based on pivot */
				// for (k=1;k<=m;++k){
					// A[i][k] = A[i][k] - A[j][k]*scale;
				// }
				// b[i] = b[i] - b[j]*scale; /* same for b */
			// }	
		// }			
	// }

	// void back_sub(double **A, double *x, double *b, int m){
	  // int i,j;
	  // x[m] = b[m]/A[m][m];

	  // for (i=m-1;i>=1;--i){
		// x[i] = b[i];
		// for (j=i+1;j<=m;++j){
		  // x[i] -= A[i][j]*x[j];
		// }
		// x[i]/=A[i][i];
	  // }
	// }
// };