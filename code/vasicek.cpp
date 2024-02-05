#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include "kalman_filter.h"
#include "vasicek.h"

/*using std::cout;
using std::endl;
using std::vector;
//using std::string;
using std::stringstream;
using std::ifstream;
using std::ios;*/
using namespace std;

vasicek::vasicek(string filename = "data/bond_data_aug2011.txt",
                 size_t total_bonds = 9, size_t total_factors = 3,
                 size_t total_prices_quoted = 5) :
    nbond(total_bonds), nfac(total_factors), nmeas(total_prices_quoted)
{
    // Allocate memory
    perform_allocations();

    // Initialize quantities
    initialize_model_parameters();
    read_data_initialize_arrays(filename);
};

vasicek::~vasicek(){};

size_t vasicek::get_nbond() { return nbond; };
size_t vasicek::get_nfac() { return nfac; };
size_t vasicek::get_nmeas() { return nmeas; };
vector<double> vasicek::get_beta() { return beta; };
double vasicek::get_alpha() { return alpha; };

void vasicek::set_beta(double &tau)
{
    for (size_t i=0; i<get_nfac(); i++)
        beta[i] = ( 1 - exp(-kappa[i]*tau) ) / kappa[i];
};

void vasicek::set_alpha(double &tau)
{
    // Formula for zero covariances
    alpha = 0;
    vector<double> bb = get_beta();
    for (size_t i=0; i<gamma.size(); i++){
        alpha += gamma[i] * (bb[i] - tau) / (kappa[i]*kappa[i]) -
                 ( sigma[i]*sigma[i] * bb[i]*bb[i] ) / (4*kappa[i]);
    };
};

void vasicek::perform_allocations()
{
    // Allocate arrays, rows only
    time_tau.resize(get_nbond());
    bond_px_last.resize(get_nbond());
    bond_px_bid.resize(get_nbond());
    bond_px_low.resize(get_nbond());
    bond_px_high.resize(get_nbond());
    bond_px_ask.resize(get_nbond());
    bond_rate_last.resize(get_nbond());
    bond_rate_bid.resize(get_nbond());
    bond_rate_low.resize(get_nbond());
    bond_rate_high.resize(get_nbond());
    bond_rate_ask.resize(get_nbond());

    zero_rate.resize(get_nbond());
    h5_time.resize(get_nbond());
    h5_px_last.resize(get_nbond());
    h5_px_bid.resize(get_nbond());
    h5_px_low.resize(get_nbond());
    h5_px_high.resize(get_nbond());
    h5_px_ask.resize(get_nbond());
    h5_bond_rate_last.resize(get_nbond());
    h5_bond_rate_bid.resize(get_nbond());
    h5_bond_rate_low.resize(get_nbond());
    h5_bond_rate_high.resize(get_nbond());
    h5_bond_rate_ask.resize(get_nbond());

    // Three-parameter Vasicek model
    beta.resize(get_nfac());

    state_model_A.resize(get_nfac());
    for (size_t i = 0; i < get_nfac(); i++)
        state_model_A[i].resize(get_nfac());

    state_error_F.resize(get_nfac());
    for (size_t i = 0; i < get_nfac(); i++)
        state_error_F[i].resize(get_nfac());

    state_covariance_Rv.resize(get_nfac());
    for (size_t i = 0; i < get_nfac(); i++)
        state_covariance_Rv[i].resize(get_nfac());

    Kalman_covariance.resize(get_nfac());
    for (size_t i = 0; i < get_nfac(); i++)
        Kalman_covariance[i].resize(get_nfac());

    state_const.resize(get_nfac());
    state_var_intrinscic.resize(get_nfac());

    // Five price quotation for measured rate
    measure_covariance_Rw.resize(get_nmeas());
    for (size_t i=0; i<get_nmeas(); i++)
        measure_covariance_Rw[i].resize(get_nmeas());

    state_measure_C.resize(get_nmeas());
    for (size_t i=0; i<get_nmeas(); i++)
        state_measure_C[i].resize(get_nfac());

    measure_const.resize(get_nmeas());
    state_var_measured.resize(get_nmeas());
};

void vasicek::set_state_eqn_matrices(double &tau)
{
    for (size_t i=0; i<get_nfac(); i++){

        state_covariance_Rv[i][i] = 0.5*sigma[i]*sigma[i]/kappa[i]*
            ( 1 - exp(-2 * kappa[i] * tau) );

        state_const[i] = theta[i] * kappa[i] * beta[i];
        state_model_A[i][i] = exp( -kappa[i] * tau );
        state_error_F[i][i] = 1.0;
    };
};

void vasicek::set_measure_eqn_matrices(size_t &ibond, size_t &step, double &tau)
{
    double tau_inv = 1.0 / tau;

    for (size_t i=0; i<get_nmeas(); i++)
        measure_const[i] = -alpha * tau_inv;

    for (size_t i=0; i<get_nmeas(); i++)
        for (size_t j=0; j<get_nfac(); j++)
            state_measure_C[i][j] = beta[j] * tau_inv;

    state_var_measured[0] = bond_rate_last[ibond][step] - measure_const[0];
    state_var_measured[1] = bond_rate_bid[ibond][step]  - measure_const[1];
    state_var_measured[2] = bond_rate_low[ibond][step]  - measure_const[2];
    state_var_measured[3] = bond_rate_high[ibond][step] - measure_const[3];
    state_var_measured[4] = bond_rate_ask[ibond][step]  - measure_const[4];

    /* Consider the measure covariance matrix diagonal
       with variances equal to 10% the measured bond rate.*/
    measure_covariance_Rw[0][0] = 0.1*bond_rate_last[ibond][step];
    measure_covariance_Rw[1][1] = 0.1*bond_rate_bid[ibond][step];
    measure_covariance_Rw[2][2] = 0.1*bond_rate_low[ibond][step];
    measure_covariance_Rw[3][3] = 0.1*bond_rate_high[ibond][step];
    measure_covariance_Rw[4][4] = 0.1*bond_rate_ask[ibond][step];
};

void vasicek::initialize_model_parameters()
{
    /* Three factor Vasicek model from
       Affine Term-Structure Models: Theory
       and Implementation, D. J. Bolder,
       Bank of Canada working paper 2001-15. */

    kappa.resize(get_nfac());
    theta.resize(get_nfac());
    sigma.resize(get_nfac());
    lambda.resize(get_nfac());
    gamma.resize(get_nfac());

    kappa[0] = 0.06; kappa[1] = 0.30 ; kappa[2] = 0.70;
    theta[0] = 0.01; theta[1] = 0.02 ; theta[2] = 0.04;
    sigma[0] = 0.02; sigma[1] = 0.05 ; sigma[2] = 0.03;
    lambda[0] = -0.20; lambda[1] = -0.50 ; lambda[2] = -0.15;
    for (size_t i=0; i<gamma.size(); i++)
        gamma[i] = kappa[i]*kappa[i]*
            (theta[i] - sigma[i]*lambda[i]/kappa[i]) - 0.5*sigma[i]*sigma[i];
};

void vasicek::read_data_initialize_arrays(string filename)
{
    size_t nrec, counter = 0;
    string temp_str1, temp_str2;
    double T_maturity, time, tau_inv;
    double px_last, px_bid, px_low, px_high, px_ask;

    ifstream input (filename.c_str(), ios::in);
    if (!input.is_open())
        {
            cout << "Could not find file " << filename << endl;
            abort();
        };

    while (!input.eof()) {
        input >> temp_str1 >> temp_str2 >> nrec >> T_maturity;

        // Update the vector containers
        if (!input.fail()) {
            bond_maturity_T.push_back(T_maturity);

            // Column creation of 2D arrays
            time_tau[counter].resize(nrec);
            bond_px_last[counter].resize(nrec);
            bond_px_bid[counter].resize(nrec);
            bond_px_low[counter].resize(nrec);
            bond_px_high[counter].resize(nrec);
            bond_px_ask[counter].resize(nrec);
            bond_rate_last[counter].resize(nrec);
            bond_rate_bid[counter].resize(nrec);
            bond_rate_low[counter].resize(nrec);
            bond_rate_high[counter].resize(nrec);
            bond_rate_ask[counter].resize(nrec);

            for (size_t i=0; i<nrec; i++){
                input >> time >> px_last >> px_bid
                      >> px_low >> px_high >> px_ask;

                time_tau[counter][i] = T_maturity - time;
                tau_inv = 1/time_tau[counter][i];

                /* Normalize the bond prices to unity
                   (quoted normalized to $100.00). */
                bond_px_last[counter][i] = px_last * 0.01;
                bond_px_bid[counter][i]  = px_bid * 0.01;
                bond_px_low[counter][i]  = px_low * 0.01;
                bond_px_high[counter][i] = px_high * 0.01;
                bond_px_ask[counter][i]  = px_ask * 0.01;
                bond_rate_last[counter][i] = - log(bond_px_last[counter][i])*tau_inv;
                bond_rate_bid[counter][i]  = - log(bond_px_bid[counter][i])*tau_inv;
                bond_rate_low[counter][i]  = - log(bond_px_low[counter][i])*tau_inv;
                bond_rate_high[counter][i] = - log(bond_px_high[counter][i])*tau_inv;
                bond_rate_ask[counter][i]  = - log(bond_px_ask[counter][i])*tau_inv;

            };

            // Bond counter (should go up to 9)
            counter++;

        };
    };

    input.close();
};

void vasicek::evolve_single_bond(size_t ibond, size_t time_steps)
{
    bool convergence;
    double zero_temp;

    // Instantiate a three-parameter Kalman filter
    kalman_filter K_filter(get_nfac(),get_nmeas());

    /* Estimate the initial Kalman covariance equal to half the identity 
       matrix. Most probably it is a bad estimate, but the filtering
       algorithm is supposed to converge it asymptotically. */
    for (size_t i=0; i<get_nfac(); i++)
        Kalman_covariance[i][i] = 0.5;

    /* Estimate the initial state variable equal to one of the measured
       rates. Remember to subtract the constant vector for each timestep
       AFTER the Kalman filtering process for that timestep. The formula
       connecting the parameters to the zero-rate is zero-rate = sum(params). */
    for (size_t i=0; i<get_nfac(); i++)
        state_var_intrinscic[i] = 0.2 * state_var_measured[i];

    // Time evolution loop
    for (size_t step=0; step<time_steps; step++){

        /* Evaluate the coefficients alpha and beta of the
           equation: bond price = exp(alpha - sum(beta_i*xbar_i))*/
        set_beta(time_tau[ibond][step]);
        set_alpha(time_tau[ibond][step]);

        // Set the matrices to the state equation
        set_state_eqn_matrices(time_tau[ibond][step]);

        // Set the matrices to the measurement equation
        set_measure_eqn_matrices(ibond, step, time_tau[ibond][step]);

        // Perform Kalman filtering
        K_filter.Kalman_driver
            (state_model_A, state_error_F, state_measure_C, state_covariance_Rv,
             measure_covariance_Rw, Kalman_covariance, state_var_measured,
             state_var_intrinscic, convergence);

        // Rescale the state variable with the constant term
        for (size_t i=0; i<get_nfac(); i++)
          state_var_intrinscic[i] -= state_const[i];
        //state_var_intrinscic[i] -= state_const[i];
        
        // Test for convergence and record Kalman estimate of zero-rate
        if (convergence) {
            zero_temp = 0;
            for (size_t i=0; i<get_nfac(); i++)
                zero_temp -= state_var_intrinscic[i];
              //zero_temp += state_var_intrinscic[i];

            zero_rate[ibond].push_back(zero_temp);

            h5_time[ibond].push_back( bond_maturity_T[ibond] - time_tau[ibond][step] );
            h5_px_last[ibond].push_back(bond_px_last[ibond][step]);
            h5_px_bid[ibond].push_back(bond_px_bid[ibond][step]);
            h5_px_low[ibond].push_back(bond_px_low[ibond][step]);
            h5_px_high[ibond].push_back(bond_px_high[ibond][step]);
            h5_px_ask[ibond].push_back(bond_px_ask[ibond][step]);
            h5_bond_rate_last[ibond].push_back(bond_rate_last[ibond][step]);
            h5_bond_rate_bid[ibond].push_back(bond_rate_bid[ibond][step]);
            h5_bond_rate_low[ibond].push_back(bond_rate_low[ibond][step]);
            h5_bond_rate_high[ibond].push_back(bond_rate_high[ibond][step]);
            h5_bond_rate_ask[ibond].push_back(bond_rate_ask[ibond][step]);
        };
            
    };
};

void vasicek::write_output()
{
    string filename;
    stringstream int2str;
    ofstream output;

    filename = "Bond_";

    for(size_t ibond=0; ibond<get_nbond(); ibond++) {
        int2str << ibond + 1;
        if (ibond != 1) {
            filename = "Bond_" + int2str.str() + ".txt";
            int2str.str("");
            output.open(filename.c_str());

            output << bond_maturity_T[ibond] << endl;

            for (size_t i=0; i<zero_rate[ibond].size(); i++)
                output << h5_time[ibond][i]           << "\t"
                       << h5_px_last[ibond][i]        << "\t"
                       << h5_px_bid[ibond][i]         << "\t"
                       << h5_px_low[ibond][i]         << "\t"
                       << h5_px_high[ibond][i]        << "\t"
                       << h5_px_ask[ibond][i]         << "\t"
                       << h5_bond_rate_last[ibond][i] << "\t"
                       << h5_bond_rate_bid[ibond][i]  << "\t"
                       << h5_bond_rate_low[ibond][i]  << "\t"
                       << h5_bond_rate_high[ibond][i] << "\t"
                       << h5_bond_rate_ask[ibond][i]  << "\t"
                       << zero_rate[ibond][i]         << endl;

            output.close();
        };
    };

    for (size_t i=0; i<10; i++) {
        int2str << i + 1;
        filename = "zero_curve_" + int2str.str() + ".txt";
        int2str.str("");
        output.open(filename.c_str());

        for(size_t ibond=0; ibond<get_nbond(); ibond++) {
            if (ibond != 1) {
                output << bond_maturity_T[ibond] << "\t"
                       << zero_rate[ibond][i]    << endl;
            };
        };

        output.close();
    };

    filename = "zero_surface.txt";
    output.open(filename.c_str());
    for(size_t ibond=0; ibond<get_nbond(); ibond++)
        if (ibond != 1) {
            for (size_t i=0; i<zero_rate[ibond].size(); i++)
                output << h5_time[ibond][i]      << "\t"
                       << bond_maturity_T[ibond] << "\t"
                       << zero_rate[ibond][i]    << endl;
        };

    output.close();
};


void vasicek::vasicek_driver()
{
    /* Loop over different bonds to record term structure.
       Skip bond 1, possible problem with data */
    evolve_single_bond(0, time_tau[0].size());
    for (size_t bonds=2; bonds<get_nbond(); bonds++){
      evolve_single_bond(bonds, time_tau[bonds].size());
    };

    // Write output to ascii format
    write_output();
};

int main()
{
    vasicek start;
    start.vasicek_driver();
};
