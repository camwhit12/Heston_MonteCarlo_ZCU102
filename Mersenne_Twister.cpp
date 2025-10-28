#include "Mersenne_Twister.h"

uint32_t mersenne_twister(mt_state* state, uint32_t seed)
{
    uint32_t* state_array = &(state->state_array[0]);
    
    // Initialize on first call
    if (state->initialized == 0) {
        state_array[0] = seed;
        
        //Grow and sow MT Array 
        for(int i=1; i<MT_n; i++)
        {
            seed = MT_f*(seed^(seed >> (MT_w-2))) + i;
            state_array[i] = seed;
        }
        
        state->state_index = 0;
        state->initialized = 1;
    }
    
    //Point to current state location 
    int k = state->state_index;

    //j = (k+1) % MT_n    
    int j = k-(MT_n-1);
    if (j<0) j += MT_n;

    //Compute x
    uint32_t x = (state_array[k] & MT_UMASK) | (state_array[j] & MT_LMASK);

    //Compute xA
    uint32_t xA = x >> 1;

    //if(x % 2) != 0
    if (x & 0x00000001UL) xA ^= MT_a;

    //Computes (k+MT_m) % MT_n
    j = k-(MT_n-MT_m);
    if (j < 0) j += MT_n;

    //Compute Twisted Value
    x = state_array[j] ^ xA;
    state_array[k++] = x;

    //Resets k to 0 when it exceeds 623
    if (k >= MT_n) k=0;
    state->state_index = k;

    //tempering
    uint32_t y = x^ (x >> MT_u);
             y = y ^ ((y << MT_s) & MT_b);
             y = y ^ ((y << MT_t) & MT_c);
    uint32_t z = y ^ (y >> MT_l);

    return z;
}
