#ifndef VASICEK_H
#define VASICEK_H

#ifndef ARRAY2D
#define ARRAY2D
using std::vector;
typedef vector<vector<double> > array2D;
#endif

using std::string;

class vasicek
{
private:
    /* nbond = number of different bonds
       nmeas = number of price quotations per bond
       nfac  = number of Vasicek factors to the rate
       tau   = time to maturity - current time (yr) */
    size_t nbond, nfac, nmeas;
    vector<double> bond_maturity_T;
    array2D time_tau;
    array2D bond_px_last, bond_px_bid, bond_px_low;
    array2D bond_px_high, bond_px_ask;
    array2D bond_rate_last, bond_rate_bid, bond_rate_low;
    array2D bond_rate_high, bond_rate_ask;

    // Model parameters
    vector<double> kappa, theta, sigma, lambda, gamma;

    /* Linear model arrays
       state_var_intrinscic = xhat (Kalman estimated)
       state_var_measured   = ybar (quotation types, e.g. last, bid...) */
    double alpha;
    array2D state_model_A, state_error_F, state_measure_C;
    array2D state_covariance_Rv, measure_covariance_Rw, Kalman_covariance;
    vector<double> state_const, measure_const, beta;
    vector<double> state_var_intrinscic, state_var_measured;

    // Quantities to write in hdf5 output
    array2D h5_px_last, h5_px_bid, h5_px_low;
    array2D h5_px_high, h5_px_ask;
    array2D h5_bond_rate_last, h5_bond_rate_bid, h5_bond_rate_low;
    array2D h5_bond_rate_high, h5_bond_rate_ask;
    array2D h5_time, zero_rate;

    // Utility functions
    void perform_allocations();
    void set_state_eqn_matrices(double &);
    void set_measure_eqn_matrices(size_t &,size_t &,double &);

public:
    // Constructor and destructor
    vasicek(string,size_t,size_t,size_t);
    ~vasicek();

    size_t get_nbond();
    size_t get_nfac();
    size_t get_nmeas();
    vector<double> get_beta();
    double get_alpha();
    void set_beta(double &);
    void set_alpha(double &);

    void initialize_model_parameters();
    void read_data_initialize_arrays(string);
    void evolve_single_bond(size_t,size_t);
    void write_output();

    void vasicek_driver();
};

#endif
