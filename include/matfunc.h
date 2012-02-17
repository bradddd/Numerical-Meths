#ifndef _matfunc_h_included_
#define _matfunc_h_included_

class matfunc{
  public:
	static void upper_triangulate(double **A, double *b, int m, int pivotFlag);
	static void back_sub(double **A, double *x, double *b, int m);
}