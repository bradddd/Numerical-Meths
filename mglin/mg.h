
#define NPRE  3
#define NPOST 3
#define NGMAX 15

void addint(double ***uf, double ***uc, double ***res, int nf);
void copy(double ***aout, double ***ain, int n);
void fill0(double ***u, int n);
void interp(double ***uf, double ***uc, int nf);
void relax(double ***u, double ***rhs, int n, double C);
void resid(double ***res, double ***u, double ***rhs, int n, double C);
void rstrct(double **uc, double **uf, int nc);
void rstrct_3d(double ***uc, double ***uf, int nc);
void slvsml(double ***u, double ***rhs, double C);
void slvsml_3d(double ***u, double ***rhs, double C);
