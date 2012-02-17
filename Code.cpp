// Brad Matthiesen
// matthiesen@uchicago.edu
// Winter 2012

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <cstring>
#include <sstream>
#include "nrutil2.h"

using namespace std;
using std::string;

// This is the basic class to contain the basic environment details
// Notably, it contains the constants observed in solving the PDE
// as well as the data member that corresponds to the temperature
// values in a 3D cartesian space.

class MatrixT{
  public:
	int nx,ny,nz;
	double l_x,l_y,l_z;
	double dx,dy,dz;
	double C_x,C_y,C_z;
	double*** data;       // Temperature data
	
	double dt, alpha;

	// Default construct makes a space of 30x30x30
	MatrixT() : nx(30), ny(30), nz(30), l_x(1), l_y(1), l_z(1) {
		dx=l_x/nx;
		dy=l_y/ny;
		dz=l_z/nz;
		allocateData();
		InitializeGauss();
		setDs();
	}

	// Constructor taking dimensions as 3 double inputs
	MatrixT(double _nx, double _ny, double _nz) : nx(_nx), ny(_ny), nz(_nz), l_x(1), l_y(1), l_z(1)  {
		dx=l_x/nx;
		dy=l_y/ny;
		dz=l_z/nz;
		allocateData();
		InitializeGauss();
		setDs();
	} 

	~MatrixT(){}

	// This calculates the d's and C's for solving the PDE
	void setDs(){
		dx=l_x/nx;
		dy=l_y/ny;
		dz=l_z/nz;
		
		alpha = .01;
		dt = .0005;

		C_x = alpha*dt/(dx*dx);
		C_y = alpha*dt/(dy*dy);
		C_z = alpha*dt/(dz*dz);
		//sdfsd
		return;
	}

	// allocate the 3D Matrix to hold T using nrutil.c function
	void allocateData(){
		data = d3tensor(1,nx,1,ny,1,nz);
	}

	void deAllocateData(){
		cout << "trying to delete" << endl;
		free_d3tensor(data,1,nx,1,ny,1,nz);
	}

	void InitializeGauss(){
					
		for ( int x = 1 ; x <= nx ; x++){
			for ( int y = 1 ; y <= ny ; y++){
				for ( int z = 1 ; z <= nz ; z++){
					//Set intial gaussian temp
					data[x][y][z] = 5*exp((-1.0)*pow((5.0*(double)x/nx)-2.5,2))
								    *exp((-1.0)*pow((5.0*(double)y/ny)-2.5,2))
								    *exp((-1.0)*pow((5.0*(double)z/nz)-2.5,2));
					//cout << data[x][y][z] << " ";
				}
			}
		}
		BoundaryCondition();	//set the boundary condition;
		return;
	}
		
	//Set BoundaryCondition To Zeros
	void BoundaryCondition(){
		/*
		// for ( int x = 0 ; x <= nx ; x++){
			// for ( int y = 0 ; y <= ny ; y++){
				// data[x][y][0] = 0;
				// data[x][y][0] = 0;
				// data[x][y][nz] = 0;
				// data[x][y][nz] = 0;
			// }
		// }
		
		// for ( int x = 0 ; x <= nx ; x++){
			// for ( int z = 0 ; z <= nz ; z++){
				// data[x][0][z] = 0;
				// data[x][0][z] = 0;
				// data[x][ny][z] = 0;
				// data[x][ny][z] = 0;
			// }
		// }
		
		// for ( int z = 0 ; z <= nz ; z++){
			// for ( int y = 0 ; y <= ny ; y++){
				// data[0][y][z] = 0;
				// data[0][y][z] = 0;
				// data[nx][y][z] = 0;
				// data[nx][y][z] = 0;
			// }
		// }
		*/
		
		for ( int x = 1 ; x <= nx ; x++){
			for ( int y = 1 ; y <= ny ; y++){
				data[x][y][1] = 0;
				data[x][y][1] = 0;
				data[x][y][nz] = 0;
				data[x][y][nz] = 0;
			}
		}
		
		for ( int x = 1 ; x <= nx ; x++){
			for ( int z = 1 ; z <= nz ; z++){
				data[x][1][z] = 0;
				data[x][1][z] = 0;
				data[x][ny][z] = 0;
				data[x][ny][z] = 0;
			}
		}
		
		for ( int z = 1 ; z <= nz ; z++){
			for ( int y = 1 ; y <= ny ; y++){
				data[1][y][z] = 0;
				data[1][y][z] = 0;
				data[nx][y][z] = 0;
				data[nx][y][z] = 0;
			}
		}
		
		return;
	}
	
	//Set BoundaryConditions be the same
	void PeriodicCondition(){
				
		for ( int x = 1 ; x <= nx ; x++){
			for ( int y = 1 ; y <= ny ; y++){
				data[x][y][1] = 0;
				data[x][y][1] = 0;
				data[x][y][nz] = 0;
				data[x][y][nz] = 0;
			}
		}
		
		for ( int x = 1 ; x <= nx ; x++){
			for ( int z = 1 ; z <= nz ; z++){
				data[x][1][z] = 0;
				data[x][1][z] = 0;
				data[x][ny][z] = 0;
				data[x][ny][z] = 0;
			}
		}
		
		for ( int z = 1 ; z <= nz ; z++){
			for ( int y = 1 ; y <= ny ; y++){
				data[1][y][z] = 0;
				data[1][y][z] = 0;
				data[nx][y][z] = 0;
				data[nx][y][z] = 0;
			}
		}
		
		return;
	}
	
	// To incorporate constant conditions
	void DirchletCondition(double (*func)(int,int,int)){
		for ( int x = 1 ; x <= nx ; x++){
			for ( int y = 1 ; y <= ny ; y++){
				for ( int z = 1 ; z <= nz ; z++){
					//Set intial gaussian temp
					data[x][y][z] += func(x,y,z);
					//cout << data[x][y][z] << " ";
				}
			}
		}
		return;
	}

	void exportT(ofstream& myFile, int slice, int time) {
		int x = slice;
		for ( int y = 0 ; y < ny ; y++) {
			for ( int z = 0 ; z < nz ; z++){
				// Export infomation from Temp Matrix for a give z=
				myFile << time << ","<< data[x][y][z] << "," << x << "," << y << "," << z << endl;
				//myFile << time << "\t"<< (*data)[x][y][z] << "\t" << x << "\t" << y << "\t" << z << endl;
			}
		}
		return;
	}
	
	//void printToFile(int p, int q, int r, double ***M) {
	void printToFile(ofstream* f) {
		//int p = nx+1, q = ny+1, r = nz+1;
		int p = nx, q = ny, r = nz;
		int prod = p*q*r;
		//int index = (int)nx/2;
		
		//cout << "prod is: " << prod << endl;
				
		
		
		//file.write((char *)(*data[index]), sizeof(double)*prod);
		f->write((char *)(&data[1][1][1]), sizeof(double)*prod);		
	}

};

// Basic functions to perform a gaussian elimination
class gausselim{
  public:
	static void upper_triangulate(double **A, double *b, int m, int pivotFlag){
	  //cout << "In upper triangular" << endl;
	  
	  int i,j,k;
	  double scale;
		
		for (j=1;j<m;++j){           /* loop over columns */
			
			// set pivot index
			i = j+1;
			
			//cout << "j is " << j << " and i is " << i << endl;
			
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

	static void back_sub(double **A, double *x, double *b, int m){
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
	
	static void mprint(double **matrix, int m, string label){
		int i, j;
		cout << label;

		for (i = 1; i <= m; ++i){
			for (j = 1; j <= m; ++j){
				printf("%10.2f ", matrix[i][j]);
			}
			printf("\n");
		}
		printf("\n------------------------\n");
	}
	static void mprint(double **matrix, int m, string label,ofstream& file){
		int i, j;
		file << label;

		for (i = 1; i <= m; ++i){
			for (j = 1; j <= m; ++j){
				file << " " << matrix[i][j];
			}
			file << endl;
		}
		file << "END MATRIX" << endl;	
	}
	static void vprint(double *vector, int m, string label){
		int i;
		cout << label;

		for (i = 1; i <= m; ++i){
			printf("%10.2f ", vector[i]);
		}

		printf("\n------------------------\n");
	}
	static void vprint(double *vector, int m, string label, ofstream& file) {
		int i;
		file << label;

		for (i = 1; i <= m; ++i){
			file << " " << vector[i];
		}

		file << "END VECTOR" << endl;
	}
};
	
class FTCS_Solver{
	public:
	
		/*
		*  General Solution for FTCS
		*
		*  T(n+1) = T ( 1 - 2*C_x - 2*C_y - 2*C_z ) 
		*		  + C_x(left and right)
		*		  + C_y(top and bottom) 
		*		  + C_z(forward and backward)
		*/
		
		string probsize;

		MatrixT* temp;
		MatrixT* new_temp;
		
		double alpha, dt;
		double dx,dy,dz;
		double C_x,C_y,C_z;
						
		FTCS_Solver(MatrixT* input) {
									
			temp = input;
			return;
		}	
		
		// no real setup like Crank, just manage overhead for exporting
		void setup(){
			dx=temp->l_x/temp->nx;
			dy=temp->l_y/temp->ny;
			dz=temp->l_z/temp->nz;

			C_x = temp->alpha*temp->dt/(dx*dx);
			C_y = temp->alpha*temp->dt/(dy*dy);
			C_z = temp->alpha*temp->dt/(dz*dz);
			
			alpha = temp->alpha;
			dt = temp->dt;
			
			new_temp = new MatrixT(temp->nx,temp->ny,temp->nz);
			cout << "declared new_temp" << endl;
			
			return;
		}
		
		// iterate through
		void solve_nextT() {		
			
			//cout << temp->data[temp->nx/2][temp->ny/2][temp->nz/2] << C_x << " " << C_y << " " << C_z << endl;
			
			for ( int x = 2 ; x < temp->nx ; x++){
				for ( int y = 2 ; y < temp->ny ; y++){
					for ( int z = 2 ; z < temp->nz ; z++){
						
						//cout << "in loop where x,y,z" << x << "," << y << "," << z << endl;    // use sparingly as it will cause putty to crash
						
						new_temp->data[x][y][z] = temp->data[x][y][z] * ( 1 - 2*C_x - 2*C_y - 2*C_z )
												 + C_x*(temp->data[x-1][y][z] + temp->data[x+1][y][z]) 
									 + C_y*(temp->data[x][y+1][z] + temp->data[x][y-1][z]) 
												 + C_z*(temp->data[x][y][z+1] + temp->data[x][y][z-1]);
						
						
					}
				}
			}
			
			MatrixT* tmp = temp;
			temp = new_temp;
			new_temp = tmp;
			
			return;			
		}
		
		void solve_nextT( double (*func)(int,int,int)) {		
			
			//cout << temp->data[temp->nx/2][temp->ny/2][temp->nz/2] << C_x << " " << C_y << " " << C_z << endl;
			
			for ( int x = 2 ; x < temp->nx ; x++){
				for ( int y = 2 ; y < temp->ny ; y++){
					for ( int z = 2 ; z < temp->nz ; z++){
						
						//cout << "in loop where x,y,z" << x << "," << y << "," << z << endl;    // use sparingly as it will cause putty to crash
						
						new_temp->data[x][y][z] = temp->data[x][y][z] * ( 1 - 2*C_x - 2*C_y - 2*C_z )
												 + C_x*(temp->data[x-1][y][z] + temp->data[x+1][y][z]) 
									 + C_y*(temp->data[x][y+1][z] + temp->data[x][y-1][z]) 
												 + C_z*(temp->data[x][y][z+1] + temp->data[x][y][z-1]) + func(x,y,z);
						
						
					}
				}
			}
			
			MatrixT* tmp = temp;
			temp = new_temp;
			new_temp = tmp;
			
			return;			
		}
				
		void solve_nextT( int val) {		
			
			//cout << temp->data[temp->nx/2][temp->ny/2][temp->nz/2] << C_x << " " << C_y << " " << C_z << endl;
			
			for ( int x = 2 ; x < temp->nx ; x++){
				for ( int y = 2 ; y < temp->ny ; y++){
					for ( int z = 2 ; z < temp->nz ; z++){
						
						//cout << "in loop where x,y,z" << x << "," << y << "," << z << endl;    // use sparingly as it will cause putty to crash
						
						new_temp->data[x][y][z] = temp->data[x][y][z] * ( 1 - 2*C_x - 2*C_y - 2*C_z )
												 + C_x*(temp->data[x-1][y][z] + temp->data[x+1][y][z]) 
									 + C_y*(temp->data[x][y+1][z] + temp->data[x][y-1][z]) 
												 + C_z*(temp->data[x][y][z+1] + temp->data[x][y][z-1]) + val;
						
						
					}
				}
			}
			
			MatrixT* tmp = temp;
			temp = new_temp;
			new_temp = tmp;
			
			return;			
		}
		
		
		double solve(){
			time_t begin, end; 
			time(&begin);
			
			std::stringstream out;
			
			out << temp->nx << "x" << temp->ny << "x" << temp->nz ;
			probsize = out.str();
			
			// populate C's and other constants
			setup();
			
			string filename = string("output_FTCS_") + probsize + string(".bin") ;
			ofstream file(filename.c_str( ), ios::binary);
									
			int ctr = 0, ctr2 = 0;
			
			for ( long t = 0 ; t <= 10000 ; t++) {
			//for ( long t = 0 ; t < 10 ; t++) {	
				//calc next t
				
				
				solve_nextT();
				//cout << "after solveNextT" << endl;
				if ( ctr2++ == 250 ) {
					//cout << "in loop where t equals " << t << " and count equals " << ctr++ << endl;
					temp->printToFile(&file);
					ctr2=1;
					ctr++;
				}
			}
			cout << "Counter equals: " << ctr << endl;
			
			
			
			file.close();
			
			time(&end);
			return difftime(end, begin);
			
		}
		
		// take a given T and export it to a file
		// deprecated to use Export function of MatrixT class

		void exportT(ofstream& myFile,int slice, int time) {
			temp->exportT(myFile,slice, time);
		}
		
};

void FTCS_Driver() {
	cout << "Matrix Solver" << endl;

	int matsize[] = {10,20,30,40,50,60,70,80};
	
	ofstream outputfile;
	outputfile.open ("FTCS_Benchmark.txt");
	  
	for ( int i = 0 ; i < 6 ; i++) {
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		cout << "Solving a system of size: " << matsize[i] << endl;
		FTCS_Solver ft1(&m1);
  		long long time = ft1.solve();
		
		outputfile << ft1.probsize << " : " << time << "secs" << endl;
	}

	outputfile.close();
	return;
}

class Crank_Solver{
	public:
		/*
		* General Solution for FTCS
		*	
		* Tn+1(i,j,k) = T(i,j,k) + (1/2)*(( C_x*( Tn+1(i+1,j,k) + Tn+1(i-1,j,k))
		*								+C_y*( Tn+1(i,j+1,k) + Tn+1(i,j-1,k))
		*								+C_z*( Tn+1(i,j,k+1) + Tn+1(i,j,k-1))
		*								-6*C_xyz*Tn+1(i,j,k)
		*								)
		*								+
		*								( C_x*( Tn(i+1,j,k) + T(i-1,j,k))
		*								+C_y*( Tn(i,j+1,k) + T(i,j-1,k))
		*								+C_z*( Tn(i,j,k+1) + Tn(i,j,k-1))
		*								-6*C_xyz*Tn(i,j,k)
		*								)
		*							)
		*
		* Tn+1(i,j,k)-(1/2)*( C_x*( Tn+1(i+1,j,k) + Tn+1(i-1,j,k))
		*				+C_y*( Tn+1(i,j+1,k) + Tn+1(i,j-1,k))
		*				+C_z*( Tn+1(i,j,k+1) + Tn+1(i,j,k-1))
		*				-6*C_xyz*Tn+1(i,j,k)
		*				)
		*= T(i,j,k) + (1/2)* ( C_x*( Tn(i+1,j,k) + T(i-1,j,k))
		*						+C_y*( Tn(i,j+1,k) + T(i,j-1,k))
		*						+C_z*( Tn(i,j,k+1) + Tn(i,j,k-1))
		*						-6*C_xyz*Tn(i,j,k)
		*						)						
		*
		*/
		string probsize; 
		
		double **A, **A_stor, *b, *x;

		double C;
		long n;
		
		double alpha, dt;
		double dx,dy,dz;
		double nx,ny,nz;
		double C_x,C_y,C_z;
		
		int b_t;
		
	  	MatrixT* temp;
		MatrixT* b_mat;
		
		Crank_Solver(MatrixT* input){
			n = (long)(input->nx) *(input->ny) *(input->nz); 
			A      = dmatrix(1,n,1,n);
			A_stor = dmatrix(1,n,1,n);
			b      = dvector(1,n);
			
			temp = input;
			b_t = 0;
			
			// populate C's and other constants
			setup();
			return;
		}	
		
		// no real setup like Crank, just manage overhead for exporting
		void setup(){
		
			nx=temp->nx; ny=temp->ny; nz=temp->nz;
			
			dx=temp->l_x/temp->nx; dy=temp->l_y/temp->ny; dz=temp->l_z/temp->nz;

			C_x = temp->alpha*temp->dt/(dx*dx);
			C_y = temp->alpha*temp->dt/(dy*dy);
			C_z = temp->alpha*temp->dt/(dz*dz);
			
			alpha = temp->alpha;
			dt = temp->dt;
			
			b_mat = new MatrixT(temp->nx,temp->ny,temp->nz);
			
			//initialize A
			initializeA();
			return;
		}
		
		void initializeA(){
			long C1,C2,C3,C4,C5,C6,C7;
			
			memset(&A[1][1],0,n*n*sizeof(double));
			
			for ( int i = 1; i < n ; i++) {
				
				C1=i-nx*ny;
				C2=i-nx;
				C3=i-1;
				C4=i;
				C5=i+1;
				C6=i+nx;
				C7=i+nx*ny;
								
				// main case
				if ( i > nx*ny && i < (n-nx*ny)) {
					//set all
					A[i][C1] = -1.0*C_z/2;
					A[i][C2] = -1.0*C_y/2;
					A[i][C3] = -1.0*C_x/2;
					A[i][C4] = 1.0+2.0*C_x+2.0*C_y+2.0*C_z;
					A[i][C5] = -1.0*C_x/2;
					A[i][C6] = -1.0*C_y/2;
					A[i][C7] = -1.0*C_z/2;
				}
				else{
					if(C1>0 && C1<n) A[i][C1] = -1.0*C_z/2;
					if(C2>0 && C2<n) A[i][C2] = -1.0*C_y/2;
					if(C3>0 && C3<n) A[i][C3] = -1.0*C_x/2;
					if(C4>0 && C4<n) A[i][C4] = 1.0+2.0*C_x+2.0*C_y+2.0*C_z;
					if(C5>0 && C5<n) A[i][C5] = -1.0*C_x/2;
					if(C6>0 && C6<n) A[i][C6] = -1.0*C_y/2;
					if(C7>0 && C7<n) A[i][C7] = -1.0*C_z/2;
				}	
				
			
			}
			
			// Copy values of A into matrix to store
			memcpy( &A_stor[1][1] , &A[1][1], n*n*sizeof(double));
			
			return;
		}
		
		
		// Reset's the A matrix after a series of Gaussian Eliminations
		void resetA(){
			// mem the data from one to another...
			memcpy(&A[1][1], &A_stor[1][1], n*n*sizeof(double));
			return;
		}
		
		// loop of calcNextT and such
		// returns the time it took in seconds to solve the
		// problem size using nsteps number of steps and 
		// writing the matrix to file every exportint steps
		double solve(int nsteps,int exportint) {
			
			// for benchmarking
			time_t begin, end; 
			time(&begin);
			
			// prob dimensions for filename
			// std::stringstream out;  
			// out << temp->nx << "x" << temp->ny << "x" << temp->nz ;
			// probsize = out.str();
						
			
			// string filename = string("output_CN_") + probsize + string(".bin") ;
			// ofstream file(filename.c_str( ), ios::binary);
									
			int ctr = 0, ctr2 = 0;
			
			// print initial gaussian
			//temp->printToFile(&file);
			
			
			for ( long t = 1 ; t <= nsteps ; t++) {
				//calc next T 
				cout << "before calc new T" << endl;
				calcNextT();
				cout << "after calc new T" << endl;
				
				
				if ( ctr2++ == exportint ) {
					//b_mat->printToFile(&file);
					//temp->printToFile(&file);
					ctr2=1;
				}
				
			}
			
			cout << "Number of Matrices Exported: " << ctr << endl;
			
			//file.close();
			
			time(&end);
			return difftime(end, begin);				
		}
		
		// generate and populate new_temp from values in temp
		void calcNextT(){
			
			temp->BoundaryCondition();
			
			///* Code to store b vector in file
			std::stringstream out;  
			out << b_t++;
			string t_string = out.str();
			
			string bfilename = string("bfile_")+t_string+string(".txt");
			ofstream bfile(bfilename.c_str( ), ios::binary); 
			//*/
			
			// recalculate b
			
			double C_xyz = (C_x+C_y+C_z)/3;     // simplified average of C_x,y,z
			for ( int x = 1 ; x <= nx ; x++){
				for ( int y = 1 ; y <= ny ; y++){
					for ( int z = 1 ; z <= nz ; z++){
												
						b_mat->data[x][y][z] = temp->data[x][y][z] + (1/2)*(bhelper(x,y,z)-6*C_xyz*temp->data[x][y][z]);		 
						
						//bfile << x << " " << y << " " << z << endl;
						//bfile << b_mat->data[x][y][z] <<  endl;
					}
				}
			}
			
			b =  &(b_mat->data[1][1][1]);
			
			gausselim::vprint(b,n,"b vector",bfile);
			gausselim::mprint(A,n,"A matrix",bfile);
			
						
			x =  &(temp->data[1][1][1]);
			
			int flag = 1;                         // 1 for pivoting
			
			gausselim::vprint(x,n,"x vector before",bfile);
			
			
			cout << "About to upper tri" << endl;
			
			gausselim::upper_triangulate(A,b,n,flag);
			
			gausselim::mprint(A,n,"A matrix after",bfile);
			gausselim::vprint(x,n,"x vector before back",bfile);
			cout << "About to back sub" << endl;
			
			gausselim::back_sub(A,x,b,n);
			
			gausselim::vprint(x,n,"x vector",bfile);
			
			bfile.close();
			
			//reset A
			resetA();
			return;
			
		}
				
		void exportT( ofstream& myFile,int slice, int time) {
			temp->exportT(myFile,slice, time);
			return;
		} 
		
	
		
		double bhelper(int x, int y, int z) {
			double val = 0;
			if ((x+1)<=nx ) val += C_x*(temp->data[x+1][y][z]);
			if ((x-1)>=1 )  val += C_x*(temp->data[x-1][y][z]);
			if ((y+1)<=ny ) val += C_y*(temp->data[x][y+1][z]);
			if ((y-1)>=1 )  val += C_y*(temp->data[x][y-1][z]);
			if ((z+1)<=nz ) val += C_z*(temp->data[x][y][z+1]);
			if ((z-1)>=1 )  val += C_z*(temp->data[x][y][z-1]);	
			return val;
		}
};

void CN_Driver() {
	cout << "Matrix Solver" << endl;

	int matsize[] = {5,10,20,30,40,50};
	
	ofstream outputfile;
	outputfile.open ("CN_Benchmark.txt");
	  
	  
	for ( int i = 0 ; i < 2 ; i++) {
		cout << "About to declare a MatrixT of size: " << matsize[i] << endl;
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		cout << "Solving a system of size: " << matsize[i] << endl;
		Crank_Solver ft1(&m1);
		cout << "After declaring Crank " << endl;
  		long long time = ft1.solve(1,1);
		cout << "After solve" << endl;
		outputfile << ft1.probsize << " : " << time << "secs" << endl;
	}
	outputfile.close();
	return;
}

int main ()
{
	CN_Driver();

	return 0;
}

