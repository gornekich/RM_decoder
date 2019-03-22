#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RM.h"
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
    uint8_t *u_vector = (uint8_t *) malloc(m * sizeof(uint8_t));
    
    fwht_transform(size, L, L_float_out);
    
    for (int i = 0; i < size; i ++) {
        if (fabs(L_float_out[i]) > fabs(maximum))
        {
            maximum = L_float_out[i];
            max_index = i;
        }
    }
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

int main()
{
    uint8_t m = 8;
    uint8_t r = 2;
    int N = (int) pow(2, m);
    uint8_t N_max = 10;
    uint8_t *res;

    uint8_t **Z = lexicographycal_order(m);
    uint8_t z_0[] = {};
    coset_t coset = calc_cosets_vectors(Z[N-1], m, Z);
    
    float L_in[] = {-18.90604434, -14.82575615, 23.96521679, 15.93369138,
25.77987553, -21.42216549, -20.18026669, 22.49947065,
18.50198709, 17.85613746, -21.07063697, -20.52454129,
15.51812208, -19.60377389, -16.14607009, 20.29571287,
-21.21955566, 10.97186776, 19.5963652 , -24.80887147,
20.03709807, 27.78993603, -15.66374362, -19.27429533,
27.37327124, -17.01451696, -14.67711293, 17.46014811,
13.22261777, 16.61466933, -21.40415448, -28.92094281,
16.00095473, 11.49470996, 17.38428036, 21.87097512,
-20.2957407 , 19.43231875, -21.4277501 , 26.56521493,
20.86018949, 23.01476261, 15.60641794, 25.41509577,
21.25342454, -23.5339833 , 16.64908225, -20.21628359,
-20.23246683, 27.89428848, -15.48452756, 19.0747349 ,
23.86397399, 24.50035704, 24.24590837, 14.75663411,
-19.01848208, 10.88037717, -20.10094528, 15.65870485,
-19.81595181, -17.49874634, -28.01658409, -24.42659849,
19.89236882, 15.76010667, 15.78674386, 20.87284238,
-20.29238608, 29.65685729, -23.66281014, 17.53596296,
-16.7657423 , -11.61258063, -30.38260857, -15.56048778,
-14.70607798, 21.52427679, -23.86076647, 16.72632508,
23.50944273, -20.92758813, 22.35812883, -18.76270324,
-23.71214189, -21.11397551, -25.23829098, -20.09068789,
-20.41332631, 22.59462429, -18.60333263, 24.90878391,
-13.36292375, -26.31322214, -15.99923996, -17.2456367 ,
24.16923642, 21.2419329 , -21.01627067, -19.10310485,
-14.79737553, 22.27686561, 16.21620872, -17.10203938,
22.73637347, 25.12208899, -24.71578262, -27.61035581,
13.8430164 , -25.76991604, -22.70312227, 22.93946168,
-13.86561928, 16.04358311, 19.1704449 , -24.48680347,
16.48068925, 19.86151151, -24.4036152 , -26.45226777,
-22.6445736 , 22.55502052, 22.56682992, -23.34984623,
-25.11341085, -22.0589073 , 23.54381923, 19.69301885,
-20.47328513, 20.45948854, 21.41786709, -26.8218989 ,
-7.83284268, -29.39357095, 17.83102005, 17.85475621,
23.32125081, -16.70783575, -28.25118917, 14.89394345,
-15.71436167, -21.6306945 , 17.29852173, 23.1475334 ,
21.95718777, 26.38740206, -22.0103702 , -21.68411376,
15.97239776, -22.0798598 , -21.77023558, 24.18317046,
-27.56661691, -17.52703075, 21.41000495, 23.15458777,
17.27077511, -19.02950207, -20.65601087, 23.93400037,
-20.23347558, 24.56212722, -25.28694154, 17.98567057,
-20.55882858, -20.39181968, -19.44571342, -19.71565648,
-15.96697967, 23.56939828, -19.3871854 , 22.27758015,
21.59624508, 20.97349644, 16.9887586 , 20.73553886,
-22.12134631, -24.74084998, -27.70116807, -17.3857926 ,
-25.50834341, 14.28011925, -19.25648364, 22.53578833,
-12.3199224 , -23.2118417 , -18.44476226, -16.64451897,
23.07194263, -17.2837062 , 20.99403934, -13.73169009,
-15.14139556, 15.89148852, -11.82393016, 35.04537069,
-24.57929792, -16.51059304, -17.39263186, -18.98473788,
15.81423065, -15.59239438, 13.87373399, -18.17957695,
-6.66653898, -16.37376323, -24.0775649 , -13.47119905,
13.60913795, 19.63590731, 24.1355491 , 14.57510942,
26.8623831 , -20.40156279, 19.5568074 , -21.83933358,
-21.00587175, -19.98091242, -23.77761928, -20.30194271,
16.6366017 , -10.13238343, 25.39679956, -24.07214291,
22.64921988, -16.62398175, -16.1507992 , 14.22237447,
16.63408624, 8.60654666, -14.54267775, -16.79763393,
18.78812094, -15.62774222, -24.16670746, 27.34171863,
-18.45954886, -16.99350783, 28.142792 , 16.83073601,
12.00304944, 21.05205509, -19.70143623, -20.25550998,
22.02231719, -24.26477602, -19.06763135, 21.62787544,
16.22938185, 12.26680597, -18.63195591, -21.0329475 ,
-17.19107769, 15.8454065 , 18.29294563, -21.51152596};
    uint8_t *L_out;
    float *L = L_in;

    //print_Z(Z, m);
    //print_M(coset.cos, N/2, 2);
    //L_out = FD_R_2(L, m-1, Z);
    //print_V(L_out, N/2);

    // float *L = (float *) calloc(N, sizeof(float));
    // L[0] = 1;
    // L[6] = 1;

    res = RPA_RM(m, r, L_in, 0.1, N_max);

    for (int i = 0; i < N; i++)
        printf("%d", res[i]);
    printf("\n");

    return 0;
}
