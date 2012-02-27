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

#endif  // UTILS_H
