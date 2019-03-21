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
int to_index(uint8_t *vector,int size){
    int index = 0;
    for (int i=1;i<size+1;i++){
        index += pow(2,(size - i)) * vector[i];
    }
    return index;
}

uint8_t *sum(uint8_t *vector_1,uint8_t *vector_2,int size){
    uint8_t *vector = (uint8_t*) malloc(size * sizeof(uint8_t));
    for (int i = 0;i < size ; i ++){
        vector[i] = (vector_1[i] + vector_2[i])%2;
    }
    return vector;
}
coset_t calc_cosets_vectors(uint8_t *z_0, uint8_t m, uint8_t** Z){
    coset_t coset ;
    int size = pow(2,m);
    int idx = 0;
    int *coset_vector = (int*) calloc(m, sizeof(int));
    int N = pow(2 , m - 1);
    coset.cos = (int**) malloc( N * sizeof(int*));
    for (int i = 0; i < N; i ++ )
    {
        coset.cos[i] = (int*) malloc(2 * sizeof(int));
    }
    coset.d = (int*) malloc (N*2*sizeof(int));
    for (int i = 0; i < 2 * N; i ++){
            int b_2 = to_index(sum(z_0,Z[i],size),size);
            if (i > b_2)
            {
                coset.d[i] = coset.d[b_2];
                continue;
            }
            coset.cos[idx][0] = i;
            coset.cos[idx][1] = b_2;
            coset.d[i] = idx;
            idx += 1;
        }
    return coset;
}

uint8_t * prod(uint8_t* vector_1, uint8_t* vector_2,int  size){
    uint8_t *vector = (uint8_t*) malloc(size * sizeof(uint8_t));
    for (int i =0 ; i < size; i++ ){
        vector[i] = (vector_1[i] + vector_2[i]) % 2;
    }
    return vector;
}

int  sum_elements(uint8_t* vector_1,int  size){
    int sum_el = 0;
    for (int i =0 ; i < size; i++ ){
        sum_el += (vector_1[i]) % 2;
    }

    return sum_el;
}

float *LLR_proj(float *L, coset_t coset, int m){
    float *L_n = (float*) malloc(pow(2,m-1) * sizeof(float));
    for (int i=0;i<pow(2,m-1);i++){
        L_n[i] = log((exp(coset.cos[i][0] + coset.cos[i][1]) + 1) / (exp(coset.cos[i][0]) + exp(coset.cos[i][1]))); 
    }
    return L_n;
}
uint8_t *FD_R_2(float *L, uint8_t m, uint8_t** Z){
    int L_h = 0;
    int maximum = 0;
    uint8_t * maximum_value_h = (uint8_t *) malloc((m+1)* sizeof(uint8_t));

    int size = pow(2,m);
    for (int  i = 0; i < size ; i ++){
        L_h = 0;
        for (int j = 0; j < size; j ++){
            uint8_t sum_u_z = sum_elements(prod(Z[i],Z[j],m + 1),m+1);
            L_h += pow(-1,sum_u_z) * L[j];
        }
        if ( abs(L_h) > abs(maximum) ){
                maximum = L_h;
                for (int k =0; k < m; k++){
                    maximum_value_h[k + 1] = Z[i][k+1];
                }
               
            }
    }
    if (maximum < 0){
        maximum_value_h[0] = 1;
    }
    uint8_t* cwd = (uint8_t*) malloc( size * sizeof(uint8_t));
    for (int  i = 0; i < size ; i ++){
        cwd[i] = (maximum_value_h[0] + sum_elements(prod(Z[i], maximum_value_h, m+1),m+1))%2;
    }
    return 0;
}

/*
uint8_t *RPA_RM(uint8_t m, uint8_t r, float *L, float theta, uint8_t N);
*/


int main()
{
    int size = 5;
    uint8_t m;
    /*uint8_t *z_0 = (int*) calloc(pow(2,m), sizeof(uint8_t));
    */
        return 0;
}