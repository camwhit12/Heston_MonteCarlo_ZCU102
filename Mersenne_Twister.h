#ifndef MERSENNE_TWISTER_H_
#define MERSENNE_TWISTER_H_

#include <stdint.h>

// Use MT_ prefix to avoid macro name collisions
#define MT_n 624
#define MT_m 397
#define MT_w 32
#define MT_r 31
#define MT_UMASK (0xffffffffUL << MT_r)
#define MT_LMASK (0xffffffffUL >> (MT_w - MT_r))
#define MT_a 0x9908b0dfUL
#define MT_u 11
#define MT_s 7
#define MT_t 15
#define MT_l 18
#define MT_b 0x9d2c5680UL
#define MT_c 0xefc60000UL
#define MT_f 1812433253UL

typedef struct 
{
    uint32_t state_array[MT_n]; //MT Array to be filled and twisted
    int state_index; //index into MT Array, 0 <= state_array <= MT_n-1
    int initialized; //flag to track if state has been initialized
} mt_state;

//Generate random number (initializes on first call if needed)
uint32_t mersenne_twister(mt_state* state, uint32_t seed);

#endif