/*
 * Solvers for linear systems 
 *
 */
#include <math.h>
#include "nrutil.h"
#include "solver.h"
#include "utils.h"


void jacobi::solve(double **A, double *x, double *b, int size, double error)
{
    int max_iterations = 1;

    double sigma = 0.0;
    double *x_old = dvector(1,size);

    // start with a guess at 0
    for (int i = 1; i <= size; i++)
    {
        x[i] = 0.0;
        x_old[i] = 0.0;
    }

    // iterate until solution converges or we hit max_iterations
    while(max_iterations--)
    {
        // go thru the rows
        for(int i = 1; i <= size; i++)
        {
            sigma = 0.0;
            for(int j = 1; j <= size; j++)
            {
                if( i != j )
                {
                    sigma = sigma + (A[i][j] * x_old[j]);
                }
            }
            
            if (A[i][i] == 0)
                std::cout << "A[i][i] is zero..." << std::endl;

            x[i] = (b[i] - sigma) / A[i][i];
        }

        // setup the next iteration
        for (int i = 1; i <= size; i++)
        {
            x_old[i] = x[i];
        }

        //printf("Iteration #%d of jacobi\n", iterations - max_iterations); 
        //vprint(x,size,"x");
        /*double error = 0.0;
          if (error <= threshold)
          break;
          */
    }
}


void gauss_seidel::solve(double **A, double *x, double *b, int size, double error)
{
    double sigma = 0.0;
    int max_iterations = 100;
    double *x_old = dvector(1,size);

    // start with a guess at 0
    for (int i = 1; i <= size; i++)
    {
        x[i] = 0.0;
        x_old[i] = 0.0;
    }

    // iterate until solution converges or we hit max_iterations
    while(max_iterations--)
    {
        // go thru the rows
        for(int i = 1; i <= size; i++)
        {
            sigma = b[i];
            for(int j = 1; j <= i - 1; j++)
            {
                sigma = sigma + (A[i][j] * x[j]);
            }
            x[i] = sigma/A[i][i];
        }

        //double error = 0.0;
        //if (error <= threshold)
        //    break;
    }
}


void sor::solve(double **A, double *x, double *b, int size, double error)
{
    double omega = 1.0;
    double sigma = 0.0;
    int max_iterations = 100;
    double *x_old = dvector(1,size);

    // start with a guess at 0
    for (int i = 1; i <= size; i++)
    {
        x[i] = 0.0;
        x_old[i] = 0.0;
    }

    // iterate until solution converges or we hit max_iterations
    while(max_iterations--)
    {
        // go thru the rows
        for(int i = 1; i <= size; i++)
        {
            for(int j = i - 1; j <= size; j++)
            {
                sigma = sigma + (A[i][j] * x[j]);
            }

            for(int j = i + 1; j <= size; j++)
            {
                sigma = sigma + (A[i][j] * x_old[j]);
            }

            x[i] = x_old[i] + omega * (((b[i] - sigma) / A[i][i]) - x_old[i]);
        }

        //double error = 0.0;
        //if (error <= threshold)
        //    break;
    }
}


// Basic functions to perform a gaussian elimination
void gausselim::solve(double **A, double *x, double *b, int size, double error)
{
    upper_triangulate(A,b,size);
    back_sub(A,x,b,size);
}

void gausselim::upper_triangulate(double **A, double *b, int m, int pivotFlag)
{
    //std::cout << "In upper triangular" << std::endl;

    int i,j,k;
    double scale;

    for (j=1;j<m;++j){           /* loop over columns */

        // set pivot index
        i = j+1;

        //std::cout << "j is " << j << " and i is " << i << std::endl;

        if (pivotFlag != 0 ){           // with pivoting

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

void gausselim::back_sub(double **A, double *x, double *b, int m)
{
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



