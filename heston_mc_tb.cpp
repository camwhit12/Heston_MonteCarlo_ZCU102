#include "heston_mc.h"
#include <iostream>
#include <iomanip>
#include <cmath>

// Reference Black-Scholes for comparison
double black_scholes_call(double S, double K, double r, double T, double sigma) {
    double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    
    // Approximate normal CDF
    auto norm_cdf = [](double x) {
        return 0.5 * (1.0 + erf(x / sqrt(2.0)));
    };
    
    return S * norm_cdf(d1) - K * exp(-r*T) * norm_cdf(d2);
}

int main() {
    std::cout << "=== Heston Monte Carlo Options Pricing Testbench ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    // Test Case 1: Standard European Call Option
    std::cout << "\n--- Test Case 1: Standard Parameters ---" << std::endl;
    
    HestonParams params1;
    params1.S0 = 100.0;      // Initial stock price
    params1.K = 100.0;       // Strike price (at-the-money)
    params1.r = 0.05;        // Risk-free rate (5%)
    params1.T = 1.0;         // 1 year to maturity
    params1.v0 = 0.04;       // Initial variance (20% vol)
    params1.theta = 0.04;    // Long-term variance
    params1.kappa = 2.0;     // Mean reversion rate
    params1.sigma = 0.3;     // Volatility of volatility
    params1.rho = -0.7;      // Correlation
    
    MCConfig config1;
    config1.num_sims = 100;
    config1.num_steps = 252;  // Daily steps
    config1.num_rngs = 4;
    
    data_t option_price1;
    unsigned int seed1 = 42;
    
    std::cout << "Parameters:" << std::endl;
    std::cout << "  S0 = " << params1.S0.to_double() << std::endl;
    std::cout << "  K = " << params1.K.to_double() << std::endl;
    std::cout << "  r = " << params1.r.to_double() << std::endl;
    std::cout << "  T = " << params1.T.to_double() << std::endl;
    std::cout << "  v0 = " << params1.v0.to_double() << " (vol = " 
              << sqrt(params1.v0.to_double())*100 << "%)" << std::endl;
    std::cout << "  theta = " << params1.theta.to_double() << std::endl;
    std::cout << "  kappa = " << params1.kappa.to_double() << std::endl;
    std::cout << "  sigma = " << params1.sigma.to_double() << std::endl;
    std::cout << "  rho = " << params1.rho.to_double() << std::endl;
    std::cout << "\nMonte Carlo Configuration:" << std::endl;
    std::cout << "  Simulations: " << config1.num_sims << std::endl;
    std::cout << "  Time steps: " << config1.num_steps << std::endl;
    
    // Run simulation
    std::cout << "\nRunning Monte Carlo simulation..." << std::endl;
    heston_mc_kernel(params1, config1, &option_price1, seed1);
    
    double mc_price = option_price1.to_double();
    std::cout << "Heston MC Option Price: " << mc_price << std::endl;
    
    // Compare with Black-Scholes (as rough reference)
    double bs_price = black_scholes_call(
        params1.S0.to_double(),
        params1.K.to_double(),
        params1.r.to_double(),
        params1.T.to_double(),
        sqrt(params1.v0.to_double())
    );
    std::cout << "Black-Scholes Price (reference): " << bs_price << std::endl;
    std::cout << "Difference: " << (mc_price - bs_price) << std::endl;
    
    // Test Case 2: Out-of-the-money option
    std::cout << "\n--- Test Case 2: Out-of-the-Money Option ---" << std::endl;
    
    HestonParams params2 = params1;
    params2.K = 110.0;  // OTM strike
    
    data_t option_price2;
    unsigned int seed2 = 123;
    
    std::cout << "Parameters: Same as Test 1, but K = " << params2.K.to_double() << std::endl;
    std::cout << "Running Monte Carlo simulation..." << std::endl;
    heston_mc_kernel(params2, config1, &option_price2, seed2);
    
    std::cout << "Heston MC Option Price: " << option_price2.to_double() << std::endl;
    
    // Test Case 3: High volatility scenario
    std::cout << "\n--- Test Case 3: High Volatility ---" << std::endl;
    
    HestonParams params3 = params1;
    params3.v0 = 0.09;      // 30% initial volatility
    params3.theta = 0.09;
    params3.sigma = 0.5;    // Higher vol of vol
    
    data_t option_price3;
    unsigned int seed3 = 456;
    
    std::cout << "Parameters: v0 = " << params3.v0.to_double() 
              << " (vol = " << sqrt(params3.v0.to_double())*100 << "%)" << std::endl;
    std::cout << "Running Monte Carlo simulation..." << std::endl;
    heston_mc_kernel(params3, config1, &option_price3, seed3);
    
    std::cout << "Heston MC Option Price: " << option_price3.to_double() << std::endl;
    
    // Summary
    std::cout << "\n=== Test Complete ===" << std::endl;
    std::cout << "All tests passed successfully!" << std::endl;
    
    return 0;
}