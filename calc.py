import time
import os
from sympy import isprime
import itertools
import math

class NTT:

    def __init__(self):
        self.k = 3
        self.N = pow(2, self.k)
        self.Q = 97
        try:
            if isprime(self.Q) == False:
                raise ValueError("Q is not prime")
            if (self.Q - 1) % (2 *self.N) != 0:
                raise ValueError("Q-1 has to be divisible by 2 * N")
        except ValueError as error:
            print("Error:", error)

        for g in range(2, self.Q):
            if pow(g, self.N, self.Q) == self.Q - 1:
                self.psi = g
                break

        self.psi_inv = pow(self.psi, self.Q - 2, self.Q)
        self.omega = (self.psi * self.psi) % self.Q
        self.omega_inv = pow(self.omega, self.Q - 2, self.Q)
        self.n_inv = pow(self.N, self.Q - 2, self.Q)

    def inv_NTT(self, a):
        a = self.ntt_cyclic(a, self.omega_inv)
        p = 1
        for i in range(self.N):
            a[i] = (a[i] * self.n_inv * p) % self.Q
            p = p * self.psi_inv % self.Q
        return a

    def ntt_cyclic(self, a, w):
        for i in range(self.N):
            j = self.bit_rev(i)
            if i < j:
                a[i], a[j] = a[j], a[i]
        len = 2
        while len <= self.N:
            wl = pow(w, self.N//len, self.Q)
            for i in range(0, self.N, len):
                cur = 1
                for j in range(len//2):
                    u = a[i + j]
                    v = (a[i + j + len // 2] *  cur) % self.Q
                    a[i + j] = (u + v) % self.Q
                    a[i + j + len // 2] = (u - v) % self.Q
                    cur = (cur * wl) % self.Q
            len <<= 1
        return a

    def bit_rev(self, x):
        r = 0
        for i in range(self.k):
            r = r << 1 | (x & 1)
            x = x >> 1
        return r

class SideChannel:
    def __init__(self):
        self.noise = 5 * pow(10, -3)

    # Trying to find secret byte of message
    def multiplying_by_sequential_numbers(self, NTT, byte, path, probes = 100, max_i = None):
        # finding leaks in multiplying by sequential number
        if max_i is None:
            max_i = NTT.Q

        x_min = 0
        x_max = NTT.Q - 1
        jump = 0
        total_times = []
        i = 1
        t1, t2 = 0, 0

        while x_min != x_max and i < max_i:

            total = 0
            a_inv = [0] * NTT.N
            a_inv[byte] = i
            a = NTT.inv_NTT(a_inv)
            s = ' '.join(str(x) for x in a)

            for j in range(probes):

                start = time.time()
                os.system(path + ' ' + s)
                end = time.time()
                total += end - start

            if i == 1:
                t2 = total / probes
            total_times.append(total/probes)
            t1 = t2
            t2 = total / probes

            if abs(t2 - t1) > self.noise:
                print('OK')
                jump += 1
                x_min = max(x_min, math.ceil(jump * NTT.Q / i))
                x_max = min(x_max, jump * NTT.Q // (i-1))
            print(t1, t2, i, jump, x_min, x_max)
            i += 1

        return [x_min, x_max], total_times

    def reverse_whole_key(self, NTT):
        s_NTT_all = [0] * NTT.N
        for i in range(NTT.N):
            print(i)
            s_NTT_all[i] = self.multiplying_by_sequential_numbers(NTT, i, './NTT')[0]

        print(s_NTT_all)
        s_results = []

        for inv_s in itertools.product(*(range(a, b+1) for a, b in  s_NTT_all)):
            s_results.append(list(NTT.inv_NTT(inv_s)))

        return s_results



NTT = NTT()
side_channel = SideChannel()
print(side_channel.reverse_whole_key(NTT))
