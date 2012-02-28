#ifndef IT_SOLVER_H 
#define IT_SOLVER_H

#include "matrix.h"
#include "nrutil.h"


class IterativeSolver
{
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

        IterativeSolver(MatrixT* input) {

            temp = input;
            return;
        }	

        // no real setup like Crank, just manage overhead for exporting
        virtual void setup(){
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
        virtual void solve_nextT() {		

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

        virtual void solve_nextT( double (*func)(int,int,int)) {		

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

        virtual void solve_nextT( int val) {		

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

double zeroFunc(int i, int j, int k) { return 0.0; }

// This class contains the basic functionality and setup
// for solving the 3d heat equation in 3 dimensions using Jacobi 	
class Jacobi : public IterativeSolver
{
    public:

        /*
         *  General Solution for FTCS
         *
         *  T(n+1) = T ( 1 - 2*C_x - 2*C_y - 2*C_z ) 
         *		  + C_x(left and right)
         *		  + C_y(top and bottom) 
         *		  + C_z(forward and backward)
         */

        Jacobi(MatrixT* input) : IterativeSolver(input) 
        {
        }	

        // iterate through
        void solve_nextT() 
        {	
            solve_nextT(zeroFunc);
        }

        void solve_nextT( double (*func)(int,int,int)) 
        {
            int max_iterations = 10;

            while (max_iterations--) 
            {
                for ( int x = 2 ; x < temp->nx ; x++)
                {
                    for ( int y = 2 ; y < temp->ny ; y++)
                    {
                        for ( int z = 2 ; z < temp->nz ; z++)
                        {
                            // just aliasing for readability
                            double ***T = new_temp->data;
                            double ***T_old = temp->data;

                            double c_term = 1/(2*C_x + 2*C_y + 2*C_z + 1);
                            //std::cout << "in loop where x,y,z" << x << "," << y << "," << z << std::endl;    // use sparingly as it will cause putty to crash
                            T[x][y][z] =  c_term * 
                                (C_x * T[x-1][y][z] + C_x * T[x+1][y][z] + 
                                 C_y * T[x][y-1][z] + C_y * T[x][y+1][z] +
                                 C_z * T[x][y][z-1] + C_z * T[x][y][z+1]) + 
                                c_term * T_old[x][y][z];
                        }
                    }
                }
            }

            int halfX = temp->nx/2; int halfY = temp->ny/2; int halfZ = temp->nz/2;
            printf("Sample temp (middle point)= %.6f\n", new_temp->data[halfX][halfY][halfZ]);

            return;			
        }

        /*void solve_nextT( int val) {		

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
        }*/
};


#endif // IT_SOLVER
