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
#include "it_solver.h"
#include "utils.h"
#include "matrix.h"

extern "C"{
#include "nrutil.h"
}





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
			
			int ctr = 0, ctr2 = 0;
			
			if (benchmark){
				for ( long t = 0 ; t <= nsteps; t++) {
					solve_nextT();
				}
			}
			else {
				
				//generate probsize string for filename
				std::stringstream out;
				out << temp->nx << "x" << temp->ny << "x" << temp->nz ;
				probsize = out.str();
				
				// This outputs the binary file containing matrices at given export intervals
				// as dictated by the input paramters
				// The file also relies upon the probsize string to help create a unique filename
				std::string filename = std::string("output_FTCS_") + probsize + std::string(".bin") ;
				std::ofstream file(filename.c_str( ), std::ios::binary);
				
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
			}
				
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
  		long long time = ft1.solve(100,1,false);
		
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
		*						+C_y*( Tn+1(i,j+1,k) + Tn+1(i,j-1,k))
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
			// calculate the size of the 2D matrix used to 
			// represent the 3D discretization
			n = (long)(input->nx) *(input->ny) *(input->nz); 

			// allocate the matrices involved
			A      = dmatrix(1,n,1,n);
			A_stor = dmatrix(1,n,1,n);
			b      = dvector(1,n);
			
			temp = input;
			b_t = 0;
			
			// populate C's and other constants
			// note: setup calls initializeA()
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
			
			std::cout << "before initialize A" << std::endl;
			initializeA();
			std::cout << "after initialize A" << std::endl;
			return;
		}
		
		void initializeA(){
			long C1,C2,C3,C4,C5,C6,C7;
			
			// this allocates the block of doubles to
			// zero where A points to
			memset(&A[1][1],0,n*n*sizeof(double));
			
			std::cout << "After memset to zero" << std::endl;
			// iterate through each row
			int count = 1;
			int numThousands = 1;

			for ( int i = 1; i < n ; i++) {
				
				C1=i-nx*ny;
				C2=i-nx;
				C3=i-1;
				C4=i;
				C5=i+1;
				C6=i+nx;
				C7=i+nx*ny;


				// Debug
				if ( count++ == 5000) {
					std::cout << "Set " << ((numThousands++)*5000) << " rows of C in A" << std::endl;
					count=0;
				}


				// This block is a set of conditionals to determine
				// if each constant should be set in a give row
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

			std::cout << "After setting C's" << std::endl;  // debug
			
			// Copy values of A into another matrix to store
			memcpy( &A_stor[1][1] , &A[1][1], n*n*sizeof(double));
			
			std::cout << "After setting additional A for storage" << std::endl;   // debug
			
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
		double solve(int nsteps,int exportint, bool benchmark) {
			
			// for benchmarking
			time_t begin, end; 
			time(&begin);
												
			int ctr = 0, ctr2 = 0;
			
			if (benchmark){
				for ( long t = 1 ; t <= nsteps; t++) {
					calcNextT(benchmark);
				}
			}
			else {
			
				//generate probsize string for filename
				std::stringstream out;
				out << temp->nx << "x" << temp->ny << "x" << temp->nz ;
				probsize = out.str();
				// This outputs the binary file containing matrices at given export intervals
				// as dictated by the input paramters
				// The file also relies upon the probsize string to help create a unique filename
				std::string filename = std::string("output_CN_") + probsize + std::string(".bin");
				std::ofstream file(filename.c_str( ), std::ios::binary);
			
				for ( long t = 1 ; t <= nsteps; t++) {
				   //calc next t
					calcNextT(benchmark);
					
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
			}
					
			time(&end);
			return difftime(end, begin);				
		}
		
		// generate and populate new_temp from values in temp
		// 
		// need to implement the benchmark 
		//
		void calcNextT(bool benchmark){
			
			temp->BoundaryCondition();
			
			if (benchmark) {
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
				x =  &(temp->data[1][1][1]);
							
				// solve using the specific solver
				s->solve(A,x,b,n);
								
				//reset A
				resetA();
			}			
			else {
				// Code to store b vector in file
				// b_t is a counter variable for the filename
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

				std::cout << "About to solve" << std::endl;
				
				// solve using the specific solver
				s->solve(A,x,b,n);
				
				bfile.close();
				
				//reset A
				resetA();
			}
			return;			
		}
				
		void exportT( std::ofstream& myFile,int slice, int time) {
			temp->exportT(myFile,slice, time);
			return;
		} 
			
		// used in recalculating the b vector
		// takes the coords for a give point in the 3D space
		// and determines its neighbors influence
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

		void deAllocateData() {
			free_dmatrix(A,1,n,1,n);
			free_dmatrix(A_stor,1,n,1,n);
			free_dvector(b,1,n);
			
			temp->deAllocateData();
			b_mat->deAllocateData();
			return;			
		}
};

void CN_Driver() {
	std::cout << "Matrix Solver" << std::endl;

	int matsize[] = {5,10,20,30,40,50};
	
	std::ofstream outputfile;
	outputfile.open ("CN_Benchmark.txt");
	 
    jacobi *j = new jacobi();
	  
	for ( int i = 3 ; i <= 3; i++) {
		std::cout << "About to declare a MatrixT of size: " << matsize[i] << std::endl;
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		std::cout << "Solving a system of size: " << matsize[i] << std::endl;
		Crank_Discretization ft1(&m1);
        ft1.setSolver(j);
		std::cout << "After Crank constructor, about to solve" << std::endl;
  		long long time = ft1.solve(10,1,true);
		std::cout << "After solve" << std::endl;
		outputfile << ft1.probsize << " : " << time << "secs" << std::endl;
		ft1.deAllocateData();
	}
	outputfile.close();
	
    delete j;
    
    return;
}

void Jacobi_Driver() {
	std::cout << "Matrix Solver" << std::endl;

	int matsize[] = {5,10,20,30,40,50};
	
	std::ofstream outputfile;
	outputfile.open ("Jacobi_Benchmark.txt");
	 
	for ( int i = 0 ; i <= 3; i++) {
		std::cout << "About to declare a MatrixT of size: " << matsize[i] << std::endl;
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		std::cout << "Solving a system of size: " << matsize[i] << std::endl;
		Jacobi j(&m1);
		std::cout << "After Jacobi constructor, about to solve" << std::endl;
  		long long time = j.solve(10,1,true);
		std::cout << "After solve" << std::endl;
		outputfile << j.probsize << " : " << time << "secs" << std::endl;
	}
	outputfile.close();
	
    return;
}

void GS_Driver() {
	
	std::cout << "Matrix Solver" << std::endl;

	int matsize[] = {5,10,20,30,40,50};
	
	std::ofstream outputfile;
	outputfile.open ("GS_Benchmark.txt");
	 
	for ( int i = 0 ; i <= 3; i++) {
		std::cout << "About to declare a MatrixT of size: " << matsize[i] << std::endl;
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		std::cout << "Solving a system of size: " << matsize[i] << std::endl;
		GaussSeidel g(&m1);
		std::cout << "After GaussSeidel constructor, about to solve" << std::endl;
  		long long time = g.solve(10,1,true);
		std::cout << "After solve" << std::endl;
		outputfile << g.probsize << " : " << time << "secs" << std::endl;
	}
	outputfile.close();
	
    return;
}

void SOR_Driver() {
	std::cout << "Matrix Solver" << std::endl;

	int matsize[] = {5,10,20,30,40,50};
	
	std::ofstream outputfile;
	outputfile.open ("SOR_Benchmark.txt");
	 
	for ( int i = 0 ; i <= 3; i++) {
		std::cout << "About to declare a MatrixT of size: " << matsize[i] << std::endl;
		MatrixT m1(matsize[i],matsize[i],matsize[i]);
		std::cout << "Solving a system of size: " << matsize[i] << std::endl;
		SOR s(&m1);
		std::cout << "After SOR constructor, about to solve" << std::endl;
  		long long time = s.solve(10,1,true);
		std::cout << "After solve" << std::endl;
		outputfile << s.probsize << " : " << time << "secs" << std::endl;
	}
	outputfile.close();
	
    return;
}

int main ()
{
	Jacobi_Driver();

	GS_Driver();

	SOR_Driver();

	return 0;
}
