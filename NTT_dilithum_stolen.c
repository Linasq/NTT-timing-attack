#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <x86intrin.h>

#define N 8
#define Q 97u
#define psi 8
#define R 128
#define Q_inv 33
#define R_inv 72
#define ITERATIONS 10000

static const int32_t zetas[8] = {
    31, 3, 44, -2, -43, 24, -36, -16
};

int32_t montgomery_reduce(int64_t a) {
    int32_t t;

    t = (int32_t)a*Q_inv & 127;
    t = (a - (int64_t)t*Q) >> 7;
    return t;
}

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

void multiply_mod(int32_t a[N], int32_t b[N], int32_t c[N])
{
    for (int i = 0; i < N; i++)
    {
        c[i] = 0;
        for (int j = 0; j < a[i]; j++)
        {
            c[i] += b[i];
        }

        while (c[i] >= Q) {
            for (int j = 0; j < Q; j++)
            {
                c[i] -= 1;
            }
        }
    }
}



int main(int argc, char **argv)
{
    int32_t a[N] = {0,0,0,0,0,0,0,0};
    int32_t s[N] = {0,0,0,0,0,0,0,50};

    int32_t c[N];

    FILE *fptr;
    fptr = fopen("statistic_cycle_measurement.csv", "w");
    //fprintf(fptr, "Some text");

    double t;
    uint32_t aux_start, aux_end;
    double results[Q][ITERATIONS];
    int32_t n[Q] = {60, 83, 35, 62, 5, 12, 43, 45, 64, 73, 53, 44, 2, 29, 7, 85, 25, 94, 23, 91, 58, 47, 31, 40, 15, 66, 56, 37, 51, 88, 78, 32, 19, 76, 46, 82, 54, 63, 13, 49, 74, 50, 81, 89, 8, 67, 18, 38, 16, 80, 48, 3, 79, 14, 22, 55, 20, 24, 1, 68, 71, 72, 87, 57, 84, 70, 26, 69, 27, 59, 4, 77, 11, 96, 52, 34, 33, 86, 6, 0, 65, 28, 92, 75, 21, 9, 30, 90, 61, 41, 42, 39, 93, 17, 95, 10, 36};
    int32_t random_sum = 0;


    for (int i = 0; i < ITERATIONS; i++)
    {

        for (int j = 0; j < 97; j++)
        {

            a[N-1] = n[j];
            // fprintf(fptr, "%d,", n[j]);

            _mm_lfence();
            uint64_t start = __rdtscp(&aux_start);
            _mm_lfence();
            multiply_mod(a, s, c);
            _mm_lfence();
            uint64_t end = __rdtscp(&aux_end);
            _mm_lfence();
            results[n[j]][i] = (double)(end - start);
            // t = (double)(end - start);
            // if(i+1 == ITERATIONS)
            //     fprintf(fptr, "%f\n", t);
            // else
            //     fprintf(fptr, "%f,", t);
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

    for (int i=0; i<Q; i++)
    {
        for (int j=0; j<ITERATIONS; j++)
        {
            if (j == ITERATIONS - 1)
                fprintf(fptr, "%f\n", results[i][j]);
            else
                fprintf(fptr, "%f,", results[i][j]);
        }
    }

    fclose(fptr);
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
