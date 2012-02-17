#ifdef __cplusplus
extern "C" {
#endif

void mprint(double **matrix, int m, char *label);
void upper_triangulate(double **A, double *b, int m, int pivotFlag);
void vprint(double *vector, int m, char *label);
void back_sub(double **A, double *x, double *b, int m);

#ifdef __cplusplus
} 
#endif