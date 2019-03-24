#ifndef _RM_DECODER_H_
#define _RM_DECODER_H_

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * Structure for cosets
 */
typedef struct {
        int **cos;
        int *d;
} coset_t;

/*
 * Decoder function
 */
uint8_t *RPA_RM(uint8_t m, uint8_t r, float *L, float theta, uint8_t N_max);

#endif //_RM_DECODER_H_