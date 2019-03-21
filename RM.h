#ifndef RM_H_
#define RM_H_

typedef struct {
        int **cos;
        int *d;
} coset_t;

float *WBAWGN(float *y, float sigma);
float sigma(float SNR);
float *AWGN(int *y, float SNR);

int to_index(uint8_t *vector, int size);
uint8_t *sum(uint8_t *vector_1, uint8_t *vector_2, int size);

/*
 * Decoder functions
 */
coset_t calc_cosets_vectors(uint8_t *z_0, uint8_t m);
float *LLR_proj(float *L, coset_t coset);
uint8_t *FD_R_2(float *L, uint8_t m);
float *divide(float *L, float c, int n);

#endif //RM_H_