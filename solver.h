#ifndef SOLVER_H 
#define SOLVER_H

/*
 * Solvers for linear systems 
 *
 */
#include <math.h>
#include "nrutil.h"
#include "utils.h"


/*
 * Abstract base class that all solvers inherit from
 *
 */
class solver {
public:
    virtual void solve(double **A, double *x, double *b, int size, double error = 0.0001) = 0;
};


class jacobi : public solver
{
public: 
    
    virtual void solve(double **A, double *x, double *b, int size, double error = 0.0001);
};


class gauss_seidel : public solver
{
public: 
    
    virtual void solve(double **A, double *x, double *b, int size, double error = 0.0001);
}; 


// Basic functions to perform a gaussian elimination
class gausselim : public solver
{
  public:
      
    void solve(double **A, double *x, double *b, int size, double error = 0.0001);
  private:
    void upper_triangulate(double **A, double *b, int m, int pivotFlag = 1);
    void back_sub(double **A, double *x, double *b, int m);

};


#endif  // SOLVER_H
