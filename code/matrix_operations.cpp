#include <iostream>
#include <cstdlib>
#include <iterator>
#include <vector>
#include "nr3.h"
#include "ludcmp.h"
#include "matrix_operations.h"

using std::cout;
using std::endl;
using std::vector;

void matrix_operations::print_matrix(array2D &mat, size_t rows, size_t cols){
    for (size_t i=0; i<rows; i++) {
        cout << endl;
        for (size_t j=0; j<cols; j++)
            cout << "mat[" << i << "][" << j << "] = " << mat[i][j] << "\t";
    };
    cout << endl << endl;
};

void matrix_operations::print_vector(vector<double> &vec){
    ostream_iterator<double> out_double (cout,"\t");
    copy (vec.begin(), vec.end(), out_double);
    cout << endl << endl;
};

void matrix_operations::mat_mat_mult(array2D &mat1, array2D &mat2, array2D &mult){
    if (mat1[0].size() != mat2.size()) {
        cout << "Ill-formed matrix multiplication" << endl;
        abort();
    };

    double sum;

    for(size_t i=0; i<mat1.size(); i++){
        for(size_t j=0; j<mat2[0].size(); j++) {
            sum = 0;
            for(size_t k=0; k<mat1[0].size(); k++)
                    sum += mat1[i][k]*mat2[k][j];
            mult[i][j] = sum;
        };
    };

};

void matrix_operations::mat_vec_mult(array2D &mat, vector<double> &vec,
                                     vector<double> &prod){
    if (mat[0].size() != vec.size()) {
        cout << "Ill-formed matrix multiplication" << endl;
        abort();
    };

    double sum;
    for(size_t i=0; i<mat.size(); i++){
        sum = 0;
        for(size_t j=0; j<vec.size(); j++)
            sum += mat[i][j]*vec[j];
        prod[i] = sum;
    };
};

void matrix_operations::matrix_inverse(array2D &input, array2D &output){
    if (input.size() != input[0].size()) {
        cout << "Non-square matrix cannot be inverted." << endl;
        abort();
    };
    if (output.size() != output[0].size()) {
        cout << "Output matrix should be square matrix." << endl;
        abort();
    };
    if (input.size()    != output.size() ||
        input[0].size() != output[0].size()) {
        cout << "Output matrix should have same dimensions as input matrix." << endl;
        abort();
    };

    MatDoub inputNR3 (input.size(),input[0].size());
    MatDoub outputNR3(input.size(),input[0].size());

    for(size_t i=0; i<input.size(); i++)
        for(size_t j=0; j<input[0].size(); j++)
            inputNR3[i][j] = input[i][j];

    LUdcmp matrix (inputNR3);
    matrix.inverse(outputNR3);

    for(size_t i=0; i<input.size(); i++)
        for(size_t j=0; j<input[0].size(); j++)
            output[i][j] = outputNR3[i][j];

};

void matrix_operations::matrix_transpose(array2D &input, array2D &output){
    if (input.size() != output[0].size()) {
        cout << "Input matrix rows should equal output matrix columns." << endl;
        abort();
    };
    if (input[0].size() != output.size()) {
        cout << "Input matrix columns should equal output matrix rows." << endl;
        abort();
    };

    for (size_t i=0; i<input.size(); i++)
        for (size_t j=0; j<input[0].size(); j++)
        output[j][i] = input[i][j];
};
