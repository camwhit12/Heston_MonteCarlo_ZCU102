#include "heston_mc.h"

// Box-Muller transform using Mersenne Twister RNG
void rng_box_muller(
    mt_state* rng_state,
    data_t& z1,
    data_t& z2
) {
    #pragma HLS INLINE off
    
    // Generate two uniform random numbers using Mersenne Twister
    uint32_t u1_int = mersenne_twister(rng_state, 0);
    uint32_t u2_int = mersenne_twister(rng_state, 0);
    
    // Convert to [0,1]
    // MT generates 32-bit integers, normalize to (0,1)
    data_t u1 = data_t(u1_int) / data_t(4294967295U);
    data_t u2 = data_t(u2_int) / data_t(4294967295U);
    
    // Ensure u1 > 0 to avoid log(0)
    if (u1 < data_t(0.000001)) u1 = data_t(0.000001);
    
    data_t r = hls::sqrt(data_t(-2.0) * hls::log(u1));
    data_t theta = data_t(2.0 * 3.14159265359) * u2;
    
    z1 = r * hls::cos(theta);
    z2 = r * hls::sin(theta);
}

// Simulate a single Heston model path
data_t heston_path_simulator(
    const HestonParams params,
    const MCConfig config,
    mt_state* rng_state
) {
    #pragma HLS INLINE off
    
    data_t S = params.S0;
    data_t v = params.v0;
    data_t dt = params.T / data_t(config.num_steps);
    data_t sqrt_dt = hls::sqrt(dt);
    
    // Simulate path
    TIME_STEPS: for (unsigned int step = 0; step < config.num_steps; step++) {
        #pragma HLS PIPELINE II=3
        #pragma HLS LOOP_TRIPCOUNT min=100 max=252 avg=252
        
        data_t z1, z2;
        rng_box_muller(rng_state, z1, z2);
        
        // Correlated Brownian motions
        data_t dW_S = z1;
        data_t dW_v = params.rho * z1 + hls::sqrt(data_t(1.0) - params.rho * params.rho) * z2;
        
        // Ensure variance stays positive
        data_t v_pos = (v > data_t(0)) ? v : data_t(0);
        data_t sqrt_v = hls::sqrt(v_pos);
        
        // Euler-Maruyama discretization for stock price
        data_t dS = params.r * S * dt + sqrt_v * S * sqrt_dt * dW_S;
        S = S + dS;
        
        // Heston variance process
        data_t dv = params.kappa * (params.theta - v_pos) * dt + 
                    params.sigma * sqrt_v * sqrt_dt * dW_v;
        v = v + dv;
    }
    
    // European call option payoff
    data_t payoff = (S > params.K) ? data_t(S - params.K) : data_t(0);
    return payoff;
}

// Main Monte Carlo kernel
void heston_mc_kernel(
    const HestonParams params,
    const MCConfig config,
    data_t* option_price,
    unsigned int seed
) {
    #pragma HLS INTERFACE mode=s_axilite port=params
    #pragma HLS INTERFACE mode=s_axilite port=config
    #pragma HLS INTERFACE mode=m_axi depth=1 port=option_price offset=slave bundle=gmem0
    #pragma HLS INTERFACE mode=s_axilite port=seed
    #pragma HLS INTERFACE mode=s_axilite port=return
    
    // Initialize Mersenne Twister state
    mt_state rng_state;
    #pragma HLS ARRAY_PARTITION variable=rng_state.state_array cyclic factor=4 dim=1
    rng_state.initialized = 0;
    
    // Initialize MT with seed (happens on first call to mersenne_twister)
    mersenne_twister(&rng_state, seed);
    
    // Accumulator for option price
    acc_t sum_payoff = 0;
    
    // Monte Carlo simulation loop
    MC_SIMS: for (unsigned int sim = 0; sim < config.num_sims; sim++) {
        #pragma HLS LOOP_TRIPCOUNT min=10000 max=100000 avg=50000
        
        data_t payoff = heston_path_simulator(params, config, &rng_state);
        sum_payoff += payoff;
    }
    
    // Calculate discounted expected payoff
    data_t discount = hls::exp(-params.r * params.T);
    data_t avg_payoff = data_t(sum_payoff) / data_t(config.num_sims);
    *option_price = discount * avg_payoff;
}
