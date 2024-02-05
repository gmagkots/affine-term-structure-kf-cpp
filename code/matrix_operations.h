#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#ifndef ARRAY2D
#define ARRAY2D
using std::vector;
typedef vector<vector<double> > array2D;
#endif

class matrix_operations
{
public:
    void print_matrix(array2D &,size_t,size_t);
    void print_vector(vector<double> &);
    void mat_mat_mult(array2D &,array2D &,array2D &);
    void mat_vec_mult(array2D &, vector<double> &, vector<double> &);
    void matrix_inverse(array2D &,array2D &);
    void matrix_transpose(array2D &,array2D &);
};
    
#endif
