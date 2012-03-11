void slvsml(double **u, double **rhs)
/* 
   Solution of the model problem on the coarsest grid, where h = 1
   2 . The right-hand side is input
   in rhs[1..3][1..3] and the solution is returned in u[1..3][1..3].
*/
{
  double h=0.5;
  //fill0(u,3);
  u[2][2] = -h*h*rhs[2][2]/4.0;
}

void slvsml_3d(double ***u, double ***rhs, double C)
/* 
   Solution of the model problem on the coarsest grid, where h = 1
   2 . The right-hand side is input
   in rhs[1..3][1..3] and the solution is returned in u[1..3][1..3].
*/
{
  void fill0(double ***u, int n);    // work for 3d?
  double h=0.5;
  fill0(u,3);
  u[2][2][2] = rhs[2][2][2]/(6*C+1);
 }
