#ifndef TYPEHELPERS_H  
#define TYPEHELPERS_H
#include "arbintrinsics.h"

void mpfr_to_arbfloat(void);

/*
Since MPFR uses more precise bit-precision we could change to
typdef{...
...int user_prec_bits; // e.g. 70
} arbfloat_t;
 
Then, handle GRT as we do already (in blocks of 64), but then 
re-round the last limb to the user precision bits.
*/

#endif