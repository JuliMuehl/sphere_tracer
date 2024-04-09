#ifndef LCG_RNG_GUARD
#define LCG_RNG_GUARD
#include <inttypes.h>
struct LCGState{
  uint32_t state;
};

static inline void seed_lcg(struct LCGState* lcg,uint32_t seed){
  lcg->state = seed;
}

static inline float next_float_lcg(struct LCGState* lcg){
  float res = (lcg->state % 10000000) / 10000000.0f;
  //See Numerical Recipes Chapter 7.1, §An Even Quicker Generator
  //See also wikipedia https://en.wikipedia.org/wiki/Linear_congruential_generator#Parameters_in_common_use
  lcg->state = (lcg->state * 1664525 + 1013904223) ;
  return res;
}
#endif
