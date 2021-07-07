#ifndef	__DWAVE__types_h
#define	__DWAVE__types_h

#include <stdint.h>

typedef signed char     	spin_t;
#define RAND_T_MAX              UINT32_MAX

#define	MAX_DEGREE	8

#define RAND_A          16807
#define PRIME           15485863

#define get_random(seed, rand)      do {        \
            rand = seed;                        \
            seed += RAND_A;                     \
            rand *= PRIME;                      \
            rand ^= (rand >> 16);               \
            rand *= PRIME;                      \
            rand ^= (rand >> 16);               \
            rand *= PRIME;                      \
            rand ^= (rand >> 16);               \
            rand *= PRIME;                      \
            rand ^= (rand >> 16);               \
} while (0);


// used in expf() approximation (with Taylor series)
#define         FACT_1_0        1.0000000000000000f
#define         FACT_1_1        1.0000000000000000f
#define         FACT_1_2        0.5000000000000000f
#define         FACT_1_3        0.1666666666666667f
#define         FACT_1_4        0.0416666666666667f

#endif	// __DWAVE_types_h
