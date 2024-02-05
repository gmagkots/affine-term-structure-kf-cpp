#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#ifndef ARRAY2D
#define ARRAY2D
using std::vector;
typedef vector<vector<double> > array2D;
#endif

#include "matrix_operations.h"

class kalman_filter
{
private:
    size_t Npar, Ndat;
    bool convergence;
    array2D state_model_matrix;           // matrix A(Npar,Npar)
    array2D state_error_matrix;           // matrix F(Npar,Npar)
    array2D state_measure_matrix;         // matrix C(Ndat,Npar)
    array2D state_vector_covariance;      // matrix Rv(Npar,Npar)
    array2D measure_vector_covariance;    // matrix Rw(Ndat,Ndat)
    
    array2D Kalman_gain;                  // matrix K(Npar,Ndat)
    array2D Kalman_covariance;            // matrix P(Npar,Npar)
    vector<double> Kalman_estimator;      // vector xhat(Npar)
    array2D Ptest, Psave;                 // Used for convergence

    matrix_operations operations;

    // Uncomment for v&v test
    // friend void test_kalman();

public:
    // Constructor and destructor
    kalman_filter(size_t,size_t);
    ~kalman_filter();

    // Get array dimensions and utility operations
    size_t get_Npar();
    size_t get_Ndat();
    matrix_operations get_operations();

    // Initialize the parameter containers
    void set_state_model_matrix(array2D &);
    void set_state_error_matrix(array2D &);
    void set_state_measure_matrix(array2D &);
    void set_vector_covariances(array2D &,array2D &);

    /* Evaluate the covariance to the Kalman
       estimate and the Kalman gain */
    void evaluate_Kalman_covariance(array2D &);
    void evaluate_Kalman_estimator(vector<double> &,vector<double> &);
    bool test_convergence();
    void Kalman_driver(array2D &,array2D &,array2D &,array2D &,array2D &,
                       array2D &,vector<double> &,vector<double> &,bool &);
};

#endif
