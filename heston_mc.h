#ifndef HESTON_MC_H
#define HESTON_MC_H

#include <ap_fixed.h>
#include <hls_math.h>
#include <hls_stream.h>
#include "Mersenne_Twister.h"

typedef ap_fixed<32,16> data_t;
typedef ap_fixed<64,32> acc_t;

struct HestonParams {
    data_t S0;  //Initial stock price
    data_t K;   //Strike price
    data_t r;   //Risk free rate
    data_t T;   //Time to maturity
    data_t v0;  //Initial variance
    data_t theta;   //Long-term variance
    data_t kappa;   //Mean reversion rate
    data_t sigma;   //Volatility of volatility (Vol of vol)
    data_t rho;     //Correlation between stock and variance
};

struct MCConfig {
    unsigned int num_sims;  //Number of simulations
    unsigned int num_steps; //Time steps per path
    unsigned int num_rngs;  //Number of parallel RNG units
};

void heston_mc_kernel(
    const HestonParams params,
    const MCConfig config,
    data_t* option_price,
    unsigned int seed
);

void rng_box_muller(
    mt_state* rng_state,
    data_t& z1,
    data_t& z2
);

data_t heston_path_simulator(
    const HestonParams params,
    const MCConfig config,
    mt_state* rng_state
);

#endif