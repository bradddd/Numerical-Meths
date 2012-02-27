#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>

#include "nrutil.h"
#include "utils.h"
#include "solver.h"

#define ERROR_TOL 1E-9

bool gauss_elim(double **A, double *b, double *x, int m);
bool upper_triangulate(double **A, double *b, int m);
bool back_sub(double **A, double *x, double *b, int m);

bool checkSolution(double *x, double *sol, int size)
{
    for (int i = 1; i <= size; i++)
    {
        if (fabs(x[i] - sol[i]) > ERROR_TOL)
        {
            printf("x[i]= %.15f sol[i]=%.15f i=%d\n", x[i], sol[i], i);
            return false;
        }
    }
    return true;
}


void solPrint(bool passed, const std::string& extra)
{
    if (passed)
        std::cout << extra << " PASSED" << std::endl;
    else
        std::cout << extra << " FAILED" << std::endl;

}

int main(int argc, char **argv){
    double **A, *b, *x, *sol;
    int n = 3;
    double threshold = 0.0001;
    A = dmatrix(1,n,1,n);
    b = dvector(1,n);
    x = dvector(1,n);
    sol = dvector(1,n);

    A[1][1] = 1; A[1][2]= 1; A[1][3] = 2;
    A[2][1] = 2; A[2][2]= 4; A[2][3] = -3;
    A[3][1] = 3; A[3][2]= 6; A[3][3] = -5;

    b[1] = 9; b[2] = 1; b[3] = 0;
    sol[1] = 1; sol[2] = 2; sol[3] = 3;
    mprint(A,n,"A original");
    vprint(b,n,"b original");
   

    std::cout << "/*" << std::endl;
    std::cout << "* TESTING Gaussian Elim Solver" << std::endl;
    std::cout << "*/" << std::endl;
    {

        solver *s = new gausselim();
        s->solve(A,x,b,n);
        mprint(A,n,"A");
        vprint(b,n,"b");
        vprint(x,n,"x");
        solPrint(checkSolution(x,sol,n), "Gaussian Elim");
    } 


    std::cout << "/*" << std::endl;
    std::cout << "* TESTING Jacobi Solver" << std::endl;
    std::cout << "*/" << std::endl;
    {
        solver *s = new jacobi();
        s->solve(A,x,b,n,0.000001);
        vprint(x,n,"x");
        solPrint(checkSolution(x,sol,n), "Jacobi");
    }
 
    
    std::cout << "/*" << std::endl;
    std::cout << "* TESTING Gauss Seidel Solver" << std::endl;
    std::cout << "*/" << std::endl;
    {
        solver * s = new gauss_seidel();
        s->solve(A,x,b,n,0.000001);
        vprint(x,n,"x");
        solPrint(checkSolution(x,sol,n), "Gauss Seidel");
    }
 
    
    std::cout << "/*" << std::endl;
    std::cout << "* TESTING  SOR Solver" << std::endl;
    std::cout << "*/" << std::endl;
    {
        solver * s = new sor();
        s->solve(A,x,b,n,0.000001);
        vprint(x,n,"x");
        solPrint(checkSolution(x,sol,n), "SOR");
    }
    
}


