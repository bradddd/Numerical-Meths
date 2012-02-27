// Brad Matthiesen
// matthiesen@uchicago.edu
// Winter 2012

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <cstring>
#include <sstream>
#include "solver.h"
#include "utils.h"

extern "C"{
#include "nrutil.h"
}



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
		std::cout << "trying to delete" << std::endl;
		free_d3tensor(data,1,nx,1,ny,1,nz);
	}

	// This function initializes the data member
	// to a 3d gaussian dist scaled accordingly
	void InitializeGauss(){
					
		for ( int x = 1 ; x <= nx ; x++){
			for ( int y = 1 ; y <= ny ; y++){
				for ( int z = 1 ; z <= nz ; z++){
					//Set intial gaussian temp for each x,y,z point
					data[x][y][z] = 5*exp((-1.0)*pow((5.0*(double)x/nx)-2.5,2))
								    *exp((-1.0)*pow((5.0*(double)y/ny)-2.5,2))
								    *exp((-1.0)*pow((5.0*(double)z/nz)-2.5,2));
				}
			}
		}
		BoundaryCondition();	//set the boundary condition;
		return;
	}
		
	//Set BoundaryCondition To Zeros
	void BoundaryCondition(){
	
		// iterate through all of the edges of the 3d region
		// and set them equal to zero
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
	// Doesn't appear to be properly implemented
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
	// This function augments the values of each point
	// in the data member by a double amount determined
	// by the function passed in by pointer
	void DirchletCondition(double (*func)(int,int,int)){
		for ( int x = 1 ; x <= nx ; x++){
			for ( int y = 1 ; y <= ny ; y++){
				for ( int z = 1 ; z <= nz ; z++){
				   // augments the value at data[x][y][z]
					// by the value returned by the function
					// passed in by pointer
					data[x][y][z] += func(x,y,z);
				}
			}
		}
		return;
	}

	// This function exports a "slice" (aka z height level)
	// from the data member and exports it to a file stream
	void exportT(std::ofstream& myFile, int slice, int time) {
		int x = slice;
		for ( int y = 0 ; y < ny ; y++) {
			for ( int z = 0 ; z < nz ; z++){
				// Export infomation from Temp Matrix for a give z=
				// 1 line of code for "csv" 2nd line for "tab-delimited"
				myFile << time << ","<< data[x][y][z] << "," << x << "," << y << "," << z << std::endl;
				//myFile << time << "\t"<< (*data)[x][y][z] << "\t" << x << "\t" << y << "\t" << z << std::endl;
			}
		}
		return;
	}
	
	//void printToFile(int p, int q, int r, double ***M) {
	void printToFile(std::ofstream* f) {
		//int p = nx+1, q = ny+1, r = nz+1;
		int p = nx, q = ny, r = nz;
		int prod = p*q*r;
		//int index = (int)nx/2;
		
		//std::cout << "prod is: " << prod << std::endl;
				
		
		
		//file.write((char *)(*data[index]), sizeof(double)*prod);
		f->write((char *)(&data[1][1][1]), sizeof(double)*prod);		
	}

};



// This class contains the basic functionality and setup
// invloved for an FTCS discretization in 3 dimensions	
class FTCS_Discretization{
	public:
	
		/*
		*  General Solution for FTCS
		*
		*  T(n+1) = T ( 1 - 2*C_x - 2*C_y - 2*C_z ) 
		*		  + C_x(left and right)
		*		  + C_y(top and bottom) 
		*		  + C_z(forward and backward)
		*/
		
		std::string probsize;

		MatrixT* temp;
		MatrixT* new_temp;
		
		double alpha, dt;
		double dx,dy,dz;
		double C_x,C_y,C_z;
						
		FTCS_Discretization(MatrixT* input) {
									
			temp = input;
			return;
		}	
		
		// no real setup like Crank, just manage overhead for exporting
		void setup(){
		   // calculates d's
			dx=temp->l_x/temp->nx;
			dy=temp->l_y/temp->ny;
			dz=temp->l_z/temp->nz;

			// calculates C's
			C_x = temp->alpha*temp->dt/(dx*dx);
			C_y = temp->alpha*temp->dt/(dy*dy);
			C_z = temp->alpha*temp->dt/(dz*dz);
			
			// sets alpha and dt
			alpha = temp->alpha;
			dt = temp->dt;

			// allocates a new MatrixT the same size as the input
			// and sets the new_temp ptr to it
			new_temp = new MatrixT(temp->nx,temp->ny,temp->nz);
			
			//std::cout << "declared new_temp" << std::endl; // used to mark where in code
			
			return;
		}
		
		// iterate through
		void solve_nextT() {		
			
			//std::cout << temp->data[temp->nx/2][temp->ny/2][temp->nz/2] << C_x << " " << C_y << " " << C_z << std::endl; // output what each C is for debugging
			
			for ( int x = 2 ; x < temp->nx ; x++){
				for ( int y = 2 ; y < temp->ny ; y++){
					for ( int z = 2 ; z < temp->nz ; z++){
						
						//std::cout << "in loop where x,y,z" << x << "," << y << "," << z << std::endl;    // use sparingly as it will cause putty to crash
						
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
			
			//std::cout << temp->data[temp->nx/2][temp->ny/2][temp->nz/2] << C_x << " " << C_y << " " << C_z << std::endl;
			
			for ( int x = 2 ; x < temp->nx ; x++){
				for ( int y = 2 ; y < temp->ny ; y++){
					for ( int z = 2 ; z < temp->nz ; z++){
						
						//std::cout << "in loop where x,y,z" << x << "," << y << "," << z << std::endl;    // use sparingly as it will cause putty to crash
						
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
			
			//std::cout << temp->data[temp->nx/2][temp->ny/2][temp->nz/2] << C_x << " " << C_y << " " << C_z << std::endl;
			
			for ( int x = 2 ; x < temp->nx ; x++){
				for ( int y = 2 ; y < temp->ny ; y++){
					for ( int z = 2 ; z < temp->nz ; z++){
						
						//std::cout << "in loop where x,y,z" << x << "," << y << "," << z << std::endl;    // use sparingly as it will cause putty to crash
						
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
		
	

		// the method to solve the system
		// returns a double equal to the number of seconds
		// it took the block to execute
		//
		// call with "true" to disable file output
		//
		// nsteps = number of "rounds" of solving to run
		// exportint = the interval to output a matrix to file
		// benchmark = bool to disable extraneous features
		double solve(int nsteps, int exportint, bool benchmark){
			time_t begin, end; 
			time(&begin);

			// populate C's and other constants
			setup();
			
			if (!benchmark){
		   	//generate probsize string for filename
				std::stringstream out;
				out << temp->nx << "x" << temp->ny << "x" << temp->nz ;
				probsize = out.str();
				
				// This outputs the binary file containing matrices at given export intervals
				// as dictated by the input paramters
				// The file also relies upon the probsize string to help create a unique filename
				std::string filename = std::string("output_FTCS_") + probsize + std::string(".bin") ;
				std::ofstream file(filename.c_str( ), std::ios::binary);
			}

			int ctr = 0, ctr2 = 0;
			
			for ( long t = 0 ; t <= nsteps; t++) {
			
			   //calc next t
				solve_nextT();
				
				//std::cout << "after solveNextT" << std::endl;
				
				if ( ctr2++ == exportint) {
					//std::cout << "in loop where t equals " << t << " and count equals " << ctr++ << std::endl;
					temp->printToFile(&file);
					ctr2=1;
					ctr++;
				}
			}
			
			// this is helpful to know how many matrices to read into matlab
			std::cout << "Number of matrices exported: " << ctr << std::endl;
			
			// will want to remove for benchmarking runs
			file.close();
		
			time(&end);
			return difftime(end, begin);
			
		}
		
		// take a given T and export it to a file
		// deprecated to use Export function of MatrixT class

		void exportT(std::ofstream& myFile,int slice, int time) {
			temp->exportT(myFile,slice, time);
		}
		
};

// This class is a driver used to benchmark runs
// using FTCS discretizations instantiated via
// the FTCS_Discretization class

void FTCS_Driver() {
	std::cout << "Matrix Solver" << std::endl;

	int matsize[] = {10,20,30,40,50,60,70,80};
	
	std::ofstream outputfile;
	outputfile.open ("FTCS_Benchmark.txt");
	  
	for ( int i = 0 ; i < 6 ; i++) {
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		std::cout << "Solving a system of size: " << matsize[i] << std::endl;
		FTCS_Discretization ft1(&m1);
  		long long time = ft1.solve();
		
		outputfile << ft1.probsize << " : " << time << "secs" << std::endl;
	}

	outputfile.close();
	return;
}


// This class contains the basic functionality and setup
// invloved for Crank-Nicholson discretization in 3 dimensions
class Crank_Discretization{
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
		std::string probsize; 
		
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

        solver *s;
		
		Crank_Discretization(MatrixT* input){
			n = (long)(input->nx) *(input->ny) *(input->nz); 
			// allocate the matrices involved
			A      = dmatrix(1,n,1,n);
			A_stor = dmatrix(1,n,1,n);
			b      = dvector(1,n);
			
			temp = input;
			b_t = 0;
			
			// populate C's and other constants
			setup();
			return;
		}	
		
		// set up the constants and at the end initialize A
		void setup(){
		
			nx=temp->nx; ny=temp->ny; nz=temp->nz;
			
			dx=temp->l_x/temp->nx; dy=temp->l_y/temp->ny; dz=temp->l_z/temp->nz;

			C_x = temp->alpha*temp->dt/(dx*dx);
			C_y = temp->alpha*temp->dt/(dy*dy);
			C_z = temp->alpha*temp->dt/(dz*dz);
			
			alpha = temp->alpha;
			dt = temp->dt;
			
			// allocates a new matrix B that will also be accessed
			// as a column vector
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
	
        void setSolver(solver *solver)
        {
            s = solver;
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
						
			
			// string filename = std::string("output_CN_") + probsize + std::string(".bin") ;
			// std::ofstream file(filename.c_str( ), std::ios::binary);
									
			int ctr = 0, ctr2 = 0;
			
			// print initial gaussian
			//temp->printToFile(&file);
			
			
			for ( long t = 1 ; t <= nsteps ; t++) {
				//calc next T 
				std::cout << "before calc new T" << std::endl;
				calcNextT();
				std::cout << "after calc new T" << std::endl;
				
				
				if ( ctr2++ == exportint ) {
					//b_mat->printToFile(&file);
					//temp->printToFile(&file);
					ctr2=1;
				}
				
			}
			
			std::cout << "Number of Matrices Exported: " << ctr << std::endl;
			
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
			std::string t_string = out.str();
			
			std::string bfilename = std::string("bfile_")+t_string+std::string(".txt");
			std::ofstream bfile(bfilename.c_str( ), std::ios::binary); 
			
			// recalculate b
			double C_xyz = (C_x+C_y+C_z)/3;     // simplified average of C_x,y,z
			for ( int _x = 1 ; _x <= nx ; _x++){
				for ( int _y = 1 ; _y <= ny ; _y++){
					for ( int _z = 1 ; _z <= nz ; _z++){
												
						b_mat->data[_x][_y][_z] = temp->data[_x][_y][_z] + (1/2)*(bhelper(_x,_y,_z)-6*C_xyz*temp->data[_x][_y][_z]);		 
						//bfile << x << " " << y << " " << z << std::endl;
						//bfile << b_mat->data[x][y][z] <<  std::endl;
					}
				}
			}
			
			b =  &(b_mat->data[1][1][1]);
			
			vprint(b,n,"b vector",bfile);
			mprint(A,n,"A matrix",bfile);
						
			x =  &(temp->data[1][1][1]);
			
			vprint(x,n,"x vector before",bfile);

			std::cout << "About to upper tri" << std::endl;
			
            s->solve(A,x,b,n);
			
			bfile.close();
			
			//reset A
			resetA();
			return;
			
		}
				
		void exportT( std::ofstream& myFile,int slice, int time) {
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
	std::cout << "Matrix Solver" << std::endl;

	int matsize[] = {5,10,20,30,40,50};
	
	std::ofstream outputfile;
	outputfile.open ("CN_Benchmark.txt");
	 
    gausselim *g = new gausselim();
    jacobi *j = new jacobi();
	  
	for ( int i = 0 ; i < 2 ; i++) {
		std::cout << "About to declare a MatrixT of size: " << matsize[i] << std::endl;
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		std::cout << "Solving a system of size: " << matsize[i] << std::endl;
		Crank_Discretization ft1(&m1);
        ft1.setSolver(j);
		std::cout << "After declaring Crank " << std::endl;
  		long long time = ft1.solve(10,1);
		std::cout << "After solve" << std::endl;
		outputfile << ft1.probsize << " : " << time << "secs" << std::endl;
	}
	outputfile.close();
	
    delete g;
    delete j;
    
    return;
}

int main ()
{
	CN_Driver();

	return 0;
}

