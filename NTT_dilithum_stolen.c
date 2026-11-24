#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <x86intrin.h>

#define N 8   // docelowo 512  , jak w FALCONie
#define Q 97u // docelowo 12289, jak w FALCONie 
#define psi 8
#define R 128
#define Q_inv 33
#define R_inv 72
#define ITERATIONS 10000

static const int32_t zetas[8] = {
    31, 3, 44, -2, -43, 24, -36, -16
};


// Const time operation
int32_t montgomery_reduce(int64_t a) {
    int32_t t;

    t = (int32_t)a * Q_inv & 127;
    t = (a - (int64_t)t * Q) >> 7;
    return t;
}


// Transform to NTT
void ntt(int32_t a[N]) {
    unsigned int len, start, j, k;
    int32_t zeta, t;

    k = 0;
    for(len = 4; len > 0; len >>= 1) {
        for(start = 0; start < N; start = j + len) {
            zeta = zetas[++k];
            for(j = start; j < start + len; ++j) {
                t = montgomery_reduce((int64_t)zeta * a[j + len]);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
        }
    }
}


// Transform from NTT
void invntt_tomont(int32_t a[N]) {
    unsigned int start, len, j, k;
    int32_t t, zeta;
    const int32_t f = 11; // mont^2/256 mod Q

    k = 8;
    for(len = 1; len < N; len <<= 1) {
        for(start = 0; start < N; start = j + len) {
            zeta = -zetas[--k];
            for(j = start; j < start + len; ++j) {
                t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = montgomery_reduce((int64_t)zeta * a[j + len]);
            }
        }
    }

    for(j = 0; j < N; ++j) {
        a[j] = montgomery_reduce((int64_t)f * a[j]);
    }
}


// Perform leaky multiplication modulo N
int multiply_mod(int32_t a, int32_t s)
{
    int32_t c = 0;
    for (int j = 0; j < a; j++)
    {
        c += s;
    }

    while (c >= Q) {
        for (int j = 0; j < Q; j++)
        {
            c -= 1;
        }
    }
    return c;
}



int main(int argc, char **argv)
{
    int32_t a[N] = {0,0,0,0,0,0,0,0};
    int32_t s[N] = {2,3,4,5,6,7,8,50};
    int32_t c[N];

    FILE *fptr;
    char filename[64];

    //double t;
    uint32_t aux_start, aux_end;
    double results[Q][ITERATIONS];
    int32_t n[Q] = {60, 83, 35, 62, 5, 12, 43, 45, 64, 73, 53, 44, 2, 29, 7, 85, 25, 94, 23, 91, 58, 47, 31, 40, 15, 66, 56, 37, 51, 88, 78, 32, 19, 76, 46, 82, 54, 63, 13, 49, 74, 50, 81, 89, 8, 67, 18, 38, 16, 80, 48, 3, 79, 14, 22, 55, 20, 24, 1, 68, 71, 72, 87, 57, 84, 70, 26, 69, 27, 59, 4, 77, 11, 96, 52, 34, 33, 86, 6, 0, 65, 28, 92, 75, 21, 9, 30, 90, 61, 41, 42, 39, 93, 17, 95, 10, 36};
    int32_t random_sum = 0;


    // iteracja po kolejnych współczynnikach a i s
    for (int jj = 0; jj < N; jj++)
    {
        snprintf(filename, sizeof filename, "multiply_results/%02d-file.csv", jj);
        fptr = fopen(filename, "w");

        // ile obliczen dla kazdego z wynikow
        for (int i = 0; i < ITERATIONS; i++)
        {

            // iteracja po kazdym mozliwym wspolczynniku
            for (int j = 0; j < Q; j++)
            {

                a[jj] = n[j];

                _mm_lfence();
                uint64_t start = __rdtscp(&aux_start);
                _mm_lfence();
                c[jj] = multiply_mod(a[jj], s[jj]);
                _mm_lfence();
                uint64_t end = __rdtscp(&aux_end);
                _mm_lfence();
                results[n[j]][i] = (double)(end - start);

                int r = rand() % 20;
                if (r & 1)
                {
                    random_sum += r;
                }
                else
                {
                    random_sum -= r;
                }

            }
        }

        for (int i = 0; i < Q; i++)
        {
            for (int j = 0; j < ITERATIONS; j++)
            {
                if (j == ITERATIONS - 1)
                    fprintf(fptr, "%f\n", results[i][j]);
                else
                    fprintf(fptr, "%f,", results[i][j]);
            }
        }
        fclose(fptr);
    }

    return 0;
    // ntt(a);
    //
    // invntt_tomont(a);
    //
    // for(int j = 0; j < N; ++j) {
    //     a[j] = montgomery_reduce(a[j]);
    // }
    //
    // for (int i = 0; i < N; i++)
    // {
    //     printf("%d ", a[i]);
    // }
    // printf("\n");
}
