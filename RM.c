#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
        int **cos;
        int *d;
} coset_t;

/*
 * Private test functions
 */
float *WBAWGN(float *y, float sigma);
float sigma(float SNR);
float *AWGN(int *y, float SNR);

/*
 * Decoder functions
 */
coset_t calc_cosets_vectors(uint8_t *z_0, uint8_t m);
float *LLR_proj(float *L, coset_t *coset);
uint8_t *FD_R_2(float *L, uint8_t m);
uint8_t *RPA_RM(uint8_t m, uint8_t r, float *L, float theta, uint8_t N);



int main()
{
        return 0;
}