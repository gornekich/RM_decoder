#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RM.h"


float *divide(float *L, float c, int n)
{
    for (int i = 0; i < n; i++)
        L[i] /= c;
    return L;
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

    for (int n = 0; n < N_max; n++)
    {
        cumuLLR = (float *) calloc(N, sizeof(float));
        // 2D array for vectors in lexicographycal order
        uint8_t ** Z;
        Z = (uint8_t **) calloc(N, sizeof(uint8_t *));

        int temp_i = 0;
        for (int i = 0; i < N; i++)
        {
            Z[i] = (uint8_t *) calloc(m, sizeof(uint8_t));
            temp_i = i;
            for (int j = m - 1; j >= 0; j--)
            {
                Z[i][j] = (uint8_t) (temp_i / (int) pow(2, j));
                temp_i = temp_i % (int) pow(2, j);
            }
        }

        for (int i = 1; i < N; i++)
        {
            coset_t coset = calc_cosets_vectors(Z[i], m);
            L_est = LLR_proj(L, coset);

            if (r != 2)
            {
                y_est = RPA_RM(m-1, r-1, L_est, theta, N_max);
            }
            else
            {
                y_est = FD_R_2(L_est, m-2);
            }

            for (int k = 0; k < N; k++)
            {
                for (int l = 0; l < m; l++)
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
            cumuLLR = divide(cumuLLR, (float)(n-1), N);
            for (int k = 0; k < N; k++)
            {
                if (fabs(cumuLLR[i] - L[i]) > theta * fabs(L[i]))
                {
                    numofchange++;
                }
            }
            L = cumuLLR;
            if (numofchange == 0)
                break;
        }
    }
    L_res = (uint8_t *) calloc(N, sizeof(uint8_t));
    for (int i = 0; i < N; i++)
        L_res[i] = L[i] < 0 ? 0 : 1;
    return L_res;
}



int main()
{
    // uint8_t m = 5;
    // uint8_t r = 1;
    // float L[2] = {0.1, 0.2};
    // float theta = 0.1;
    // uint8_t N_max = 10;
    // int N = pow(2, m);

    // uint8_t *Z;
    // Z = RPA_RM(m, r, L, theta, N_max);

    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = m - 1; j >= 0; j--)
    //     {
    //         printf("%d", Z[i][);
    //     }
    //     printf("\n");
    // }    

    return 0;
}
