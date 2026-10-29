# NTT-timing-attack

Naszym problemem badawczym będzie próba odzyskania wielomianu $s(x)$ znając wielomian $a(x)$ za pomocą pomiarów czasu. Aktualnie mamy napisany kod w C, którego zadaniem jest pobranie wejścia od użytkownika, którym są współczynniki $a(x)$. Następnie program wykonuje transformację NTT obu wielomianów $a(x)$ oraz $s(x)$, uzyskując ich reprezentacje w dziedzinie NTT, a następnie przeprowadza mnożenie w tej dziedzinie. Dla testów pracujemy aktualnie z mniejszymi współczynnikami, to jest $N=8,\ \ Q=97$.

Funkcja, którą próbujemy zaatakować to `leaky_mulmod`, która wykonuje operacje mnożenia modulo w nieoptymalny sposób.
