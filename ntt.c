// Kompilacja: gcc -O1 -o ntt ntt.c
// (-O1 / -O0, żeby kompilator nie "naprawił" pętli redukcji), ale powinno juz działać bez tego

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 8            // dla wersji testowej - N=8, czyli zapisujemy na 3 bitach
#define Q 97u          // Q-1 = 2^9 * 15, więc istnieje pierwiastek 512-ego stopnia z 1 // Q = 2N * x + 1 && Q is prime

static uint32_t psi, psi_inv, omega, omega_inv, n_inv;
static uint32_t secret_hat[N];     // tajny klucz s w dziedzinie NTT

// pod sleep, dla testów
struct timespec ts = {
    .tv_sec = 0,
    .tv_nsec = 10000000  // 0.01 sekundy = 10 ms
};

/* ---------- arytmetyka modularna ---------- */

static uint32_t pow_mod(uint32_t b, uint32_t e) {
    uint64_t r = 1, x = b % Q;
    while (e) {
        if (e & 1) 
            r = r * x % Q;
        x = x * x % Q;
        e >>= 1;
    }
    return (uint32_t)r;
}


// static uint32_t pow_2(uint32_t b, uint32_t e) {
//     uint64_t r = 1;
//     while (e) { if (e & 1) r = r * b; b = b * b; e >>= 1; }
//     return (uint32_t)r;
// }

static inline uint32_t leaky_reduce(uint32_t x) {
    //printf("%u %u\n",x,Q);
    while (x >= Q) {
      //printf("a\n");
      x -= Q;
      nanosleep(&ts, NULL); // do celów testowych, nie pojawi się w implementacji
    }
    return x;
}

static inline uint32_t leaky_mulmod(uint32_t a, uint32_t b) {
    //printf("%u %u\n",a,b);
    uint32_t c = 0; // dodanie zmiennej i liczenie w ten sposób pozwala na na pominiecie argumentu -o1 przy kompilacji
    for(int i = 0; i < b; i++) {
      c += a;
    }
    return leaky_reduce(c);
}

// poprawna wersja używana w publicznych częściach NTT
static inline uint32_t mulmod(uint32_t a, uint32_t b) {
    return (uint32_t)((uint64_t)a * b % Q);
}

// static inline uint32_t mul(uint32_t a, uint32_t b) {
//     return (uint32_t)((uint64_t)a * b);
// }

/* ---------- NTT ---------- */

// odwracanie bitów - z 1010 dostajemy 0101
static uint32_t bitrev(uint32_t x, int bits) {
    uint32_t r = 0;
    for (int i = 0; i < bits; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

// cykliczne NTT rozmiaru N z pierwiastkiem w (omega)
static void ntt_cyclic(uint32_t *a, uint32_t w) {
    for (uint32_t i = 0; i < N; i++) {
        uint32_t j = bitrev(i, 3); // liczba bitów zależna od długości N
        if (i < j) { 
            uint32_t t = a[i];
            a[i] = a[j];
            a[j] = t;
        }
    }
    for (uint32_t len = 2; len <= N; len <<= 1) {
        uint32_t wl = pow_mod(w, N / len);
        for (uint32_t i = 0; i < N; i += len) {
            uint32_t cur = 1;
            for (uint32_t j = 0; j < len / 2; j++) {
                uint32_t u = a[i + j];
                uint32_t v = mulmod(a[i + j + len / 2], cur);
                a[i + j] = (u + v) % Q;
                a[i + j + len / 2] = (u + Q - v) % Q;
                cur = mulmod(cur, wl);
            }
        }
    }
}

// negacykliczne NTT: skręcenie psi^i + cykliczne NTT
static void ntt(uint32_t *a) {
    uint32_t p = 1;
    for (int i = 0; i < N; i++) { 
        a[i] = mulmod(a[i], p);
        p = mulmod(p, psi);
    }
    ntt_cyclic(a, omega);
}

// operacja odwrotna do ntt
static void intt(uint32_t *a) {
    ntt_cyclic(a, omega_inv);
    uint32_t p = 1;
    for (int i = 0; i < N; i++) {
        a[i] = mulmod(mulmod(a[i], n_inv), p);
        p = mulmod(p, psi_inv);
    }
}

/* ---------- inicjalizacja i klucz ---------- */

void init_params(void) {
    for (uint32_t g = 2; g < Q; g++) {
        //uint32_t c = pow_mod(g, 2 * N);
        if (pow_mod(g, N) == Q - 1) { psi = g; break; }   // rząd dokładnie 512 (w przyszlosci taki rzad bedzie ;)
    }
    psi_inv   = pow_mod(psi, Q - 2);
    omega     = mulmod(psi, psi);
    omega_inv = pow_mod(omega, Q - 2);
    n_inv     = pow_mod(N, Q - 2);
}


// generuje tajny klucz s o małych współczynnikach {-1,0,1}
void keygen(uint32_t seed, int16_t *s_out) {
    srand(seed);
    uint32_t s[N];
    for (int i = 0; i < N; i++) {
        int v = (rand() % 3) - 1;
        s_out[i] = (int16_t)v;
        s[i] = (v < 0) ? Q + v : v;
    }
    ntt(s);
    for (int i = 0; i < N; i++) secret_hat[i] = s[i];
}

/* ---------- podatna funkcja: r = a * s mod (x^N+1, Q) ---------- */

void poly_mul_secret(const uint16_t *a_in, uint16_t *r_out) {
    uint32_t a[N];
    
    for (int i = 0; i < N; i++) a[i] = a_in[i] % Q;
    ntt(a);
    //for (int i = 0; i < N; i++) printf("%u ", a[i]);
    //printf("\n");
    //uint32_t result;
    for (int i = 0; i < N; i++) {
        //result = leaky_mulmod(a[i], secret_hat[i]);
        //printf("%d: a[i] %u, sh[i] %u, lmul[i] %u \n", i, a[i], secret_hat[i], result);
        a[i] = leaky_mulmod(a[i], secret_hat[i]);
    }

    intt(a);
    for (int i = 0; i < N; i++) r_out[i] = (uint16_t)a[i];
}


// wypisywac do pliku?
// printowac cokolwiek?
//
int main(int argc, char **argv) {
    //int samples = (argc > 1) ? atoi(argv[1]) : 5000;

    init_params();

    int16_t s_true[N];
    keygen(1234, s_true);

    // zapisanie klucza w celu potwierdzenia ataku
    //FILE *fk = fopen("secret_true.txt", "w");
    //for (int i = 0; i < N; i++) fprintf(fk, "%d\n", s_true[i]);
    //fclose(fk);
    /* w celu szybszego wygenerowania odpowiednich zmiennych
     * powinno byc w pythonie planowo
    for(int i = 1; i < 10; i++) {
        uint32_t b[N] = {0, 0, 0, 0, 0, 0, 0, i};
        intt(b);
        for (int i = 0; i < N; i++) printf("%u ", b[i]);
        printf("\n"); 
    }
    */

    uint16_t a[N], r[N];

    if (argc != N+1) {
        printf("Niepoprawna ilość zmiennych.\nOtrzymano: %d, zamiast %d\n",argc-1,N);
        exit(1);
    }
    for (int i = 0; i < N; i++) a[i] = atoi(argv[i+1]);
    poly_mul_secret(a, r);
    /* wypisanie a oraz r. Pytanie czy printujemy, czy nie - dodatkowy szum + spam w konsoli (ale możliwy do ominięcia)
    for (int i = 0; i < N; i++) printf("%u ", a[i] % Q);
    printf("\n");
    */
    for (int i = 0; i < N; i++) printf("%u ", r[i]);
    printf("\n");
    
   
    return 0;
}

/*
85 1 8 68 27 13 90 6 
73 2 16 35 54 35 72 56 
61 3 24 2 81 57 54 9 
49 4 32 66 11 79 36 59 
37 5 40 33 38 4 18 12 
25 6 48 0 65 26 0 0 
13 7 56 64 92 48 79 15 
1 8 64 31 22 70 61 65 
86 9 72 95 49 92 43 18 
*/
