#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RM.h"
#include "RM_words.h"
#include "RM_w.h"
#include <string.h>


float *divide(float *L, float c, int n)
{
    for (int i = 0; i < n; i++)
        L[i] /= c;
    return L;
}

/*
 * Decoder functions
 */
int to_index(uint8_t *vector,int size){
    int index = 0;
    for (int i=1;i<size+1;i++)
    {
        index += pow(2,(size - i)) * vector[i-1];
    }
    return index;
}

uint8_t *sum(uint8_t *vector_1,uint8_t *vector_2,int size)
{
    uint8_t *vector = (uint8_t*) malloc(size * sizeof(uint8_t));
    for (int i = 0;i < size ; i ++)
    {
        vector[i] = (vector_1[i] + vector_2[i])%2;
    }
    return vector;
}

coset_t calc_cosets_vectors(uint8_t *z_0, uint8_t m, uint8_t** Z)
{
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
            int b_2 = to_index(sum(z_0,Z[i],m),m);
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

uint8_t * prod(uint8_t* vector_1, uint8_t* vector_2,int  size)
{
    uint8_t *vector = (uint8_t*) malloc(size * sizeof(uint8_t));
    for (int i =0 ; i < size; i++ ){
        vector[i] = (vector_1[i] * vector_2[i]) % 2;
    }
    return vector;
}

int  sum_elements(uint8_t* vector_1,int  size){
    int sum_el = 0;

    for (int i =0 ; i < size; i++ ){
        sum_el += vector_1[i];
    }
    return sum_el % 2;
}

float *LLR_proj(float *L, coset_t coset, int m)
{
    float *L_n = (float*) malloc(pow(2,m-1) * sizeof(float));
    for (int i=0;i<pow(2,m-1);i++){
        L_n[i] = log((exp(L[coset.cos[i][0]] + L[coset.cos[i][1]]) + 1) / (exp(L[coset.cos[i][0]]) + exp(L[coset.cos[i][1]]))); 
    }
    return L_n;
}

uint8_t *FD_R_2(float *L, uint8_t m, uint8_t** Z)
{
    int L_h = 0;
    int size = pow(2,m);
    uint8_t * maximum_value_h = (uint8_t *) calloc((m+1), sizeof(uint8_t));
    float *L_float_out = (float *) malloc(size*sizeof(float));
    float maximum =0;
    int max_index = 0;
    
    fwht_transform(size, L, L_float_out);
    
    for (int i = 0; i < size; i ++) {
        if (fabs(L_float_out[i]) > fabs(maximum))
        {
            maximum = L_float_out[i];
            max_index = i;
        }
    }
    free(L_float_out);
    for (int j = m - 1; j >= 0; j--)
    {
        maximum_value_h[m-j] = (uint8_t) (max_index / (int) pow(2, j));
        max_index = max_index % (int) pow(2, j);
    }

    if (maximum < 0){
        maximum_value_h[0] = 1;
    }
    uint8_t* cwd = (uint8_t*) malloc( size * sizeof(uint8_t));
    for (int  i = 0; i < size ; i ++){
        cwd[i] = (maximum_value_h[0] + sum_elements(prod(Z[i], maximum_value_h, m+1),m+1))%2;
    }
    free(maximum_value_h);
    return cwd;
}

void fwht_transform(int n, const float *src, float *dst)
{
    float adata[n];
    float bdata[n];
    float *a = adata;
    float *b = bdata;
    float *tmp;
    memcpy(a, src, sizeof(float)*n);
    
    // Fast Walsh Hadamard Transform.
    int i, j, s;
    for (i = n>>1; i > 0; i>>=1) {
        for (j = 0; j < n; j++) {
            s = j/i%2;
            b[j]=a[(s?-i:0)+j]+(s?-1:1)*a[(s?0:i)+j];
        }
        tmp = a; a = b; b = (float *)tmp;
    }
    
    memcpy(dst, a, sizeof(int)*n);
}

void print_Z(uint8_t **Z, uint8_t m)
{
    int N = pow(2,m);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf("%d", Z[i][j]);
        }
        printf("\n");
    }
    return;
}

void print_V(uint8_t *v, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%d, ", v[i]);
    }
    printf("\n");
    return;
}

void print_F(float *v, int n)
{
    for (int i = 0; i < n; ++i)
    {
        printf("%f, ", v[i]);
    }
    printf("\n");
    return;
}

void print_M(int **M, uint8_t n, uint8_t m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf("%d, ", M[i][j]);
        }
        printf("\n");
    }
    return;
}

uint8_t *RPA_RM(uint8_t m, uint8_t r, float *L, float theta, uint8_t N_max)
{
    // Code length
    int N = pow(2, m);
    // 
    float *L_est = NULL;
    //
    uint8_t *y_est = NULL;
    //
    int idx_z = 0;
    //
    float *cumuLLR = NULL;
    //
    int numofchange = 0;
    //
    uint8_t *L_res = NULL;

    uint8_t ** Z;
    Z = (uint8_t **) calloc(N, sizeof(uint8_t *));

    int temp_i = 0;
    for (int i = 0; i < N; i++)
    {
        Z[i] = (uint8_t *) calloc(m, sizeof(uint8_t));
        temp_i = i;
        for (int j = m - 1; j >= 0; j--)
        {
            Z[i][m-j-1] = (uint8_t) (temp_i / (int) pow(2, j));
            temp_i = temp_i % (int) pow(2, j);
        }
    }
    
    for (int n = 0; n < N_max; n++)
    {
        cumuLLR = (float *) calloc(N, sizeof(float));
        // 2D array for vectors in lexicographycal order

        for (int i = 1; i < N; i++)
        {
            coset_t coset = calc_cosets_vectors(Z[i], m, Z);
            L_est = LLR_proj(L, coset, m);

            if (r != 2)
            {
                y_est = RPA_RM(m-1, r-1, L_est, theta, N_max);
            }
            else
            {
                y_est = FD_R_2(L_est, m-1, Z);
            }

            for (int k = 0; k < N; k++)
            {
                idx_z = k;

                if (y_est[coset.d[idx_z]] == 0)
                {
                    cumuLLR[idx_z] += L[to_index(sum(Z[k], Z[i], m), m)];
                }
                else
                {
                    cumuLLR[idx_z] -= L[to_index(sum(Z[k], Z[i], m), m)];
                }
            }
            for (int i = 0; i < N/2; i++)
            {
                    free(coset.cos[i]);
            }
            free(coset.cos);
            free(coset.d);
        }

        numofchange = 0;
        cumuLLR = divide(cumuLLR, (float)(N-1), N);
        for (int k = 0; k < N; k++)
        {
            if (fabs(cumuLLR[k] - L[k]) > theta * fabs(L[k]))
            {
                numofchange++;
            }
        }
        L = cumuLLR;
        if (numofchange == 0)
            break;
        
    }
    L_res = (uint8_t *) calloc(N, sizeof(uint8_t));
    for (int i = 0; i < N; i++)
        L_res[i] = L[i] < 0 ? 0 : 1;

    for (int i = 0; i < N; i++)
    {
        free(Z[i]);
    }
    free(Z);
    free(cumuLLR);
    return L_res;
}

/*
 * Private functions for tests
 */
uint8_t **lexicographycal_order(uint8_t m)
{
    int N = pow(2,m);
    uint8_t ** Z;
    Z = (uint8_t **) calloc(N, sizeof(uint8_t *));

    int temp_i = 0;
    for (int i = 0; i < N; i++)
    {
        Z[i] = (uint8_t *) calloc(m, sizeof(uint8_t));
        temp_i = i;
        for (int j = m - 1; j >= 0; j--)
        {
            Z[i][m-j-1] = (uint8_t) (temp_i / (int) pow(2, j));
            temp_i = temp_i % (int) pow(2, j);
        }
    }
    return Z;
}

// /*
//  * Channel imitation functions
//  */
// float *WBAWGN(float y, float sigma)
// {
//     return -2 * y / pow(sigma, 2);
// }

// float sigma_f(float SNR)
// {
//     return sqrt(1 / pow(10, SNR/10) / 2);
// }

// float AWGN_channel(float y, float SNR)
// {
//     return -(2*y - 1); //TODO complete
//}

// def AWGN_channel(y,SNR):
//     return -(2*y - 1) + np.random.normal(0,sigma(SNR),y.shape[0])

int main()
{
     uint8_t m = 8;
     uint8_t r = 2;
     int N = (int) pow(2, m);
     uint8_t N_max = 5;
     uint8_t *res;

//     uint8_t **Z = lexicographycal_order(m);
//     uint8_t z_0[] = {};
//     coset_t coset = calc_cosets_vectors(Z[N-1], m, Z);
    

    // uint8_t *L_out;
    // float *L = L_in;

    //print_Z(Z, m);
    //print_M(coset.cos, N/2, 2);
    //L_out = FD_R_2(L, m-1, Z);
    //print_V(L_out, N/2);
    
//     float L_in[] = {-0.84291681 ,-2.90970972, 2.15106808, -1.5533445 , 1.9395258, 4.14138797,
// 1.19830696 ,-1.97100217};
//     float *L = L_in;

//     res = RPA_RM(m, r, L, 0.1, N_max);

//     for (int i = 0; i < N; i++)
//         printf("%d", res[i]);
//     printf("\n");

    int e = 0;
    for (int i = 0; i < 9991; i++)
    {
       res = RPA_RM(m, r, rm_codewords[i], 0.1, N_max);
       for (int j = 0; j < 255; ++j)
       {
           if (res[j] != rm_words[i][j])
           {
                e += 1;
                break;
            }
        }
    }

    printf("%f\n", (float)e/9991);
    return 0;
}
