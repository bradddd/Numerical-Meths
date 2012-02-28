#ifndef MATRIX_H 
#define MATRIX_H

#include "nrutil.h"
#include "util.h"

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

#endif  // MATRIX_H
