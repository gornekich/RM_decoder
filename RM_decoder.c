#include "RM_decoder.h"

/*
 * Private supplementary functions
 */

/*
 * Divide vector on constant
 */
static float *divide(float *L, float c, int n)
{
    for (int i = 0; i < n; i++)
        L[i] /= c;
    return L;
}

/*
 * Sum two vectors
 */
static uint8_t *sum(uint8_t *vector_1,uint8_t *vector_2,int size)
{
    uint8_t *vector = (uint8_t*) malloc(size * sizeof(uint8_t));
    for (int i = 0;i < size ; i ++)
    {
        vector[i] = (vector_1[i] + vector_2[i])%2;
    }
    return vector;
}

/*
 * Get index from vector
 */
static int to_index(uint8_t *vector,int size){
    int index = 0;
    for (int i=1;i<size+1;i++)
    {
        index += pow(2,(size - i)) * vector[i-1];
    }
    return index;
}

/*
 * Multiply 2 vectors
 */
static uint8_t * prod(uint8_t* vector_1, uint8_t* vector_2,int  size)
{
    uint8_t *vector = (uint8_t*) malloc(size * sizeof(uint8_t));
    for (int i =0 ; i < size; i++ ){
        vector[i] = (vector_1[i] * vector_2[i]) % 2;
    }
    return vector;
}

/*
 * Sum elements in vector
 */
static int sum_elements(uint8_t* vector_1,int  size){
    int sum_el = 0;

    for (int i =0 ; i < size; i++ ){
        sum_el += vector_1[i];
    }
    return sum_el % 2;
}

/*
 * Fast Walsh Hadamard transform
 */
static void fwht_transform(int n, const float *src, float *dst)
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

/*
 * Private decoder functions
 */

/*
 * Calculate cosets
 */
static coset_t calc_cosets_vectors(uint8_t *z_0, uint8_t m, uint8_t** Z)
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

/*
 * Calculate LLR projection
 */
static float *LLR_proj(float *L, coset_t coset, int m)
{
    float *L_n = (float*) malloc(pow(2,m-1) * sizeof(float));
    for (int i=0;i<pow(2,m-1);i++){
        L_n[i] = log((exp(L[coset.cos[i][0]] + L[coset.cos[i][1]]) + 1) / (exp(L[coset.cos[i][0]]) + exp(L[coset.cos[i][1]]))); 
    }
    return L_n;
}

/*
 * Trivial RM decoder implementation
 */
static uint8_t *FD_R_2(float *L, uint8_t m, uint8_t** Z)
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

/*
 * Reed-Muller decoder function.
 * params: m,r - RM parameters
 *         L - input vector from channel
 *         theta - threshold
 *         N_max - max iterration number
 */
uint8_t *RPA_RM(uint8_t m, uint8_t r, float *L, float theta, uint8_t N_max)
{
    // Code length
    int N = pow(2, m);
    float *L_est = NULL;
    uint8_t *y_est = NULL;
    int idx_z = 0;
    float *cumuLLR = NULL;
    int numofchange = 0;
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
            free(y_est);
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

        //free(L);
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
