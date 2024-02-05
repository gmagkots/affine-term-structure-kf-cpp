#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "matrix_operations.h"
#include "kalman_filter.h"

using std::cout;
using std::endl;
using std::vector;
using std::max;
using std::setprecision;
using std::ostream_iterator;

kalman_filter::kalman_filter(size_t paramaters_dim, size_t measure_dim) :
    Npar (paramaters_dim), Ndat (measure_dim)
{
    // Allocate arrays
    state_model_matrix.resize(get_Npar());         // rows creation
    for (size_t i = 0; i < get_Npar(); i++)
        state_model_matrix[i].resize(get_Npar());  // columns creation

    state_error_matrix.resize(get_Npar());
    for (size_t i = 0; i < get_Npar(); i++)
        state_error_matrix[i].resize(get_Npar());

    state_measure_matrix.resize(get_Ndat());
    for (size_t i = 0; i < get_Ndat(); i++)
        state_measure_matrix[i].resize(get_Npar());

    state_vector_covariance.resize(get_Npar());
    for (size_t i = 0; i < get_Npar(); i++)
        state_vector_covariance[i].resize(get_Npar());

    measure_vector_covariance.resize(get_Ndat());
    for (size_t i = 0; i < get_Ndat(); i++)
        measure_vector_covariance[i].resize(get_Ndat());

    Kalman_gain.resize(get_Npar());
    for (size_t i = 0; i < get_Npar(); i++)
        Kalman_gain[i].resize(get_Ndat());

    Kalman_covariance.resize(get_Npar());
    for (size_t i = 0; i < get_Npar(); i++)
        Kalman_covariance[i].resize(get_Npar());

    Kalman_estimator.resize(get_Npar());

    Psave.resize(get_Npar());
    for (size_t i = 0; i < get_Npar(); i++)
        Psave[i].resize(get_Npar());

    Ptest.resize(get_Npar());
    for (size_t i = 0; i < get_Npar(); i++)
        Ptest[i].resize(get_Npar());
};
kalman_filter::~kalman_filter(){};

size_t kalman_filter::get_Npar() { return Npar; };
size_t kalman_filter::get_Ndat() { return Ndat; };

matrix_operations kalman_filter::get_operations(){ return operations; };
    
void kalman_filter::set_state_model_matrix(array2D &input)
{
    for (size_t i = 0; i < get_Npar(); i++)
        for (size_t j = 0; j < get_Npar(); j++)
            state_model_matrix.at(i).at(j) = input[i][j];
};

void kalman_filter::set_state_error_matrix(array2D &input)
{
    for (size_t i = 0; i < get_Npar(); i++)
        for (size_t j = 0; j < get_Npar(); j++)
            state_error_matrix.at(i).at(j) = input[i][j];
};

void kalman_filter::set_state_measure_matrix(array2D &input)
{
    for (size_t i = 0; i < get_Ndat(); i++)
        for (size_t j = 0; j < get_Npar(); j++)
            state_measure_matrix.at(i).at(j) = input[i][j];
};

void kalman_filter::set_vector_covariances(array2D &Rv, array2D &Rw)
{
    for (size_t i = 0; i < get_Npar(); i++)
        for (size_t j = 0; j < get_Npar(); j++)
            state_vector_covariance.at(i).at(j) = Rv[i][j];

    for (size_t i = 0; i < get_Ndat(); i++)
        for (size_t j = 0; j < get_Ndat(); j++)
            measure_vector_covariance.at(i).at(j) = Rw[i][j];
};

void kalman_filter::evaluate_Kalman_covariance(array2D &covariance)
{
    array2D temp1(get_Npar(),get_Npar()),temp2(get_Npar(),get_Npar());
    array2D temp3(get_Npar(),get_Npar()),temp4(get_Npar(),get_Npar());
    array2D temp5(get_Npar(),get_Ndat()),temp6(get_Ndat(),get_Ndat());
    array2D c_tr(get_Npar(),get_Ndat()),pc_tr(get_Npar(),get_Ndat());
    array2D K_tr(get_Ndat(),get_Npar());

    // C^tr -> c_tr, P.C^tr -> pc_tr
    get_operations().matrix_transpose(state_measure_matrix,c_tr);
    get_operations().mat_mat_mult(covariance,c_tr,pc_tr);
    // C.(P.C^tr) -> temp6
    get_operations().mat_mat_mult(state_measure_matrix,pc_tr,temp6);
    // (Rw + C.P.C^tr)^-1 -> temp6
    for (size_t i = 0; i < get_Ndat(); i++)
        for (size_t j = 0; j < get_Ndat(); j++)
            temp6[i][j] = temp6[i][j] + measure_vector_covariance[i][j];
    get_operations().matrix_inverse(temp6,temp6);
    // A.(P.C^tr) -> temp5
    get_operations().mat_mat_mult(state_model_matrix,pc_tr,temp5);
    // Kalman gain K = A.P.C^tr.(Rw + C.P.C^tr)^-1
    get_operations().mat_mat_mult(temp5,temp6,Kalman_gain);

    // K.C -> temp1
    get_operations().mat_mat_mult(Kalman_gain,state_measure_matrix,temp1);
    // A - K.C -> temp2, (A - K.C)^tr -> temp1
    for (size_t i = 0; i < get_Npar(); i++)
        for (size_t j = 0; j < get_Npar(); j++)
            temp2[i][j] = state_model_matrix[i][j] - temp1[i][j];
    get_operations().matrix_transpose(temp2,temp1);
    // (A - K.C).P -> temp3
    get_operations().mat_mat_mult(temp2,covariance,temp3);
    // (A - K.C).P.(A - K.C)^tr -> temp2
    get_operations().mat_mat_mult(temp3,temp1,temp2);
    // F.Rv -> temp1
    get_operations().mat_mat_mult(state_error_matrix,state_vector_covariance,temp1);
    // F^tr -> temp3
    get_operations().matrix_transpose(state_error_matrix,temp3);
    // (F.Rv).F^tr -> temp4
    get_operations().mat_mat_mult(temp1,temp3,temp4);
    // (A - K.C).P.(A - K.C)^tr + F.Rv.F^tr -> temp1
    for (size_t i = 0; i < get_Npar(); i++)
        for (size_t j = 0; j < get_Npar(); j++)
            temp1[i][j] = temp2[i][j] + temp4[i][j];
    // K.Rw -> temp5
    get_operations().mat_mat_mult(Kalman_gain,measure_vector_covariance,temp5);
    // K^tr -> K_tr
    get_operations().matrix_transpose(Kalman_gain,K_tr);
    // (K.Rw).K^tr -> temp2
    get_operations().mat_mat_mult(temp5,K_tr,temp2);
    // Kalman covariance P = (A - K.C).P.(A - K.C)^tr + F.Rv.F^tr + K.Rw.K^tr
    for (size_t i = 0; i < get_Npar(); i++)
        for (size_t j = 0; j < get_Npar(); j++){
            Psave[i][j] = covariance[i][j];
            covariance[i][j] = temp1[i][j] + temp2[i][j];
        };
    Kalman_covariance = covariance;
};

void kalman_filter::evaluate_Kalman_estimator(vector<double> &estimator,
                                              vector<double> &measurement)
{
    vector<double> temp1(get_Npar()),temp2(get_Ndat()),temp3(get_Npar());

    // y - C.xhat -> temp2
    get_operations().mat_vec_mult(state_measure_matrix,estimator,temp2);
    for (size_t i = 0; i < get_Ndat(); i++)
        temp2[i] = measurement[i] - temp2[i];
    // K.(y - C.xhat) -> temp3
    get_operations().mat_vec_mult(Kalman_gain,temp2,temp3);
    // A.xhat -> temp1
    get_operations().mat_vec_mult(state_model_matrix,estimator,temp1);
    // Kalman estimator xhat = B + A.xhat + K.(y - C.xhat)
    for (size_t i = 0; i < get_Npar(); i++)
        estimator[i] = temp3[i] + temp1[i];
    Kalman_estimator = estimator;
};

bool kalman_filter::test_convergence()
{
    convergence = true;
    for (size_t i = 0; i < get_Npar(); i++)
        for (size_t j = 0; j < get_Npar(); j++){
            Ptest[i][j] = ( Kalman_covariance[i][j] - Psave[i][j]) /
                max( Kalman_covariance[i][j], Psave[i][j] );
            if ( abs(Ptest[i][j]) > abs(0.1*Psave[i][j]) )
                convergence = false;
        };

    return convergence;
};

void kalman_filter::Kalman_driver
(array2D &A, array2D &F, array2D &C, array2D &Rv, array2D &Rw, array2D &P,
 vector<double> &ybar, vector<double> &xhat, bool &converge)
{
    set_state_model_matrix(A);
    set_state_error_matrix(F);
    set_state_measure_matrix(C);
    set_vector_covariances(Rv,Rw);

    evaluate_Kalman_covariance(P);
    evaluate_Kalman_estimator(xhat,ybar);
    converge = test_convergence();
};



/*
void test_kalman(){
    size_t npar = 2, ndat = 2;
    kalman_filter afilter(npar,ndat);
    vector<double> xhat(npar),ybar(ndat);
    array2D A(npar,npar),F(npar,npar),C(ndat,npar),P(npar,npar);
    array2D Rv(npar,npar),Rw(ndat,ndat);
    ostream_iterator<double> out_double (cout,"\t");

    // state variable equation
    A[0][0] = 5; A[0][1] = 3; A[1][0] = 8; A[1][1] = 10;
    F[0][0] = 1; F[0][1] = 0; F[1][0] = 0; F[1][1] = 3;
    xhat[0] = 10; xhat[1] = 5;
    Rv[0][0] = 1; Rv[0][1] = 0; Rv[1][0] = 0; Rv[1][1] = 2;

    // measurement variable equation
    C[0][0] = 7; C[0][1] = 4; C[1][0] = 9; C[1][1] = 6;
    ybar[0] = 2; ybar[1] = 3;
    Rw[0][0] = 2; Rw[0][1] = 0; Rw[1][0] = 0; Rw[1][1] = 2;

    P[0][0] = 10; P[0][1] = 10; P[1][0] = 10; P[1][1] = 10;

    cout << "xhat initial is " << endl;
    copy (xhat.begin(), xhat.end(), out_double);
    cout << endl << endl;

    afilter.set_state_model_matrix(A);
    afilter.set_state_error_matrix(F);
    afilter.set_state_measure_matrix(C);
    afilter.set_vector_covariances(Rv,Rw);

    for (size_t i = 0; i < 2; i++){
        afilter.evaluate_Kalman_covariance(P);
        afilter.evaluate_Kalman_estimator(xhat,ybar);

        cout << "Kalman gain is";
        afilter.get_operations().print_matrix
            (afilter.Kalman_gain,afilter.get_Npar(),afilter.get_Ndat());

        cout << "Kalman covariance is";
        afilter.get_operations().print_matrix(P,npar,npar);

        cout << "xhat final is " << endl;
        copy (xhat.begin(), xhat.end(), out_double);
        cout << endl << endl;
    };
};

int main()
{
    test_kalman();
    return 0;
};
*/

/*The test above should output these results:

First iteration
Kalman gain K
mat[0][0] = 0.254188    mat[0][1] = 0.34662
mat[1][0] = 0.571924    mat[1][1] = 0.779896
Updated Kalman covariance P
mat[0][0] = -1.3693     mat[0][1] = -0.830928
mat[1][0] = -0.830928   mat[1][1] = -19.8696
Updated Kalman estimate xhat
xhat[0] = 2.07683 xhat[1] = -11.5771

Second iteration
Kalman gain K
mat[0][0] = 0.306052    mat[0][1] = 0.298979
mat[1][0] = -1.30572    mat[1][1] = 2.42421
Updated Kalman covariance P
mat[0][0] = 1.40571     mat[0][1] = -0.49765
mat[1][0] = -0.49765    mat[1][1] = 66.9868
Updated Kalman estimate xhat
xhat[0] = 2.06485 xhat[1] = -12.8987

The test_Kalman prototype needs to be uncommented
in kalman_filter.h*/
