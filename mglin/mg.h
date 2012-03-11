
#define NPRE  3
#define NPOST 3
#define NGMAX 15

void addint(double **uf, double **uc, double **res, int nf);
void copy(double **aout, double **ain, int n);
void fill0(double **u, int n);
void interp(double **uf, double **uc, int nf);
void relax(double **u, double **rhs, int n);
void resid(double **res, double **u, double **rhs, int n);
void rstrct(double **uc, double **uf, int nc);
void slvsml(double **u, double **rhs);
