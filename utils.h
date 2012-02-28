#ifndef UTILS_H 
#define UTILS_H

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>


inline void mprint(double **matrix, int m, std::string label){
    int i, j;
    std::cout << label;

    for (i = 1; i <= m; ++i){
        for (j = 1; j <= m; ++j){
            printf("%10.2f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n------------------------\n");
}

template<typename stream>
inline void mprint(double **matrix, int m, std::string label,stream& file){
    int i, j;
    file << label;

    for (i = 1; i <= m; ++i){
        for (j = 1; j <= m; ++j){
            file << " " << matrix[i][j];
        }
        file << std::endl;
    }
    file << "END MATRIX" << std::endl;	
}

inline void vprint(double *vector, int m, std::string label){
    int i;
    std::cout << label;

    for (i = 1; i <= m; ++i){
        printf("%10.2f ", vector[i]);
    }

    printf("\n------------------------\n");
}

template<typename stream>
inline void vprint(double *vector, int m, std::string label, stream& file) {
    int i;
    file << label;

    for (i = 1; i <= m; ++i){
        file << " " << vector[i];
    }

    file << "END VECTOR" << std::endl;
}


template<typename T>
inline void mcopy(T*** dest, T*** source, int xSize, int ySize, int zSize)
{
	for ( int i = 1 ; i < xSize; i++)
	{
		for ( int j = 1 ; j < ySize ; j++)
		{
			for ( int k = 1 ; k < zSize ; k++)
			{
				dest[i][j][k] = source[i][j][k];
			}
		}
	}
}

#endif  // UTILS_H
