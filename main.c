#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "RM_SNR.h"
#include "RM_words.h"
#include "RM_decoder.h"

/*
 * Private functions for tests
 */

/*
 * Generate matrix size 2^m x m in lexicographycal order
 */
static uint8_t **lexicographycal_order(uint8_t m)
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

/*
 * Prints byte matrix m x m
 */
static void print_Z(uint8_t **Z, uint8_t m)
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

/*
 * Prints matrix n x m
 */
static void print_M(int **M, uint8_t n, uint8_t m)
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

/*
 * Prints byte vector
 */
static void print_V(uint8_t *v, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%d, ", v[i]);
    }
    printf("\n");
    return;
}

/*
 * Prints float vector
 */
static void print_F(float *v, int n)
{
    for (int i = 0; i < n; ++i)
    {
        printf("%f, ", v[i]);
    }
    printf("\n");
    return;
}

/*
 * Main test functions
 */
static void test_rm_decoder(uint8_t m, uint8_t r, float input[9991][256],
                            uint8_t codewords[9991][256], float theta,
                            uint8_t N_max)
{
    int N = (int) pow(2, m);
    int e = 0;
    int i = 0;
    uint8_t *res;

    for (i = 0; i < 8000; i++)
    {
        printf("errors: %d, err prob:%f, i: %d \n", e, (float)e/i, i);
        res = RPA_RM(m, r, input[i], theta, N_max);
        for (int j = 0; j < N; ++j)
        {
           if (res[j] != codewords[i][j])
           {
                e += 1;
                break;
            }
        }
        if (e > 100)
            break;
    }
    free(res);

    printf("%f, i: = %d \n", (float)e/i, i);
}

int main()
{
    uint8_t m = 8;
    uint8_t r = 2;
    uint8_t N_max = 3;
    float theta = 0.1;

    /*
     * Call this function with different input vectors
     * Input vectors: rm_snr_4, rm_snr_5, rm_snr_6, rm_snr_6_5
     */
    test_rm_decoder(m, r, rm_snr_6, rm_words, theta, N_max);
    
    return 0;
}
