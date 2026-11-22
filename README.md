# NTT-timing-attack

Naszym problemem badawczym będzie próba odzyskania wielomianu $s(x)$ znając wielomian $a(x)$ za pomocą pomiarów czasu.

## Podatność i zbieranie danych

Stworzona przez nas funkcja mnożenia wielomianów w postaci NTT w C jest podatna na wyciek czasowy ze względu na sposób realizacji operacji modulo. Działanie to jest wykonywane poprzez odejmowanie $Q$ razy liczby $1$ aż współczynnik $c_i$ nie będzie należał do zbioru $\{0,1,2,3,\ldots,Q-1\}$.

W celu zebrania danych do analizy przygotowaliśmy sekretny wielomian $s$, który w postaci NTT jest reprezentowany jako $s = [0,0,0,0,0,0,0,50]$. Następnie wykonywaliśmy mnożenie go przez wielomian $a$, który w postaci NTT jest reprezentowany jako $a_i = [0,0,0,0,0,0,0,i]$, gdzie $i \in \{0,1,,2,\ldots,Q-1\}$. Dla każdego $i$ wykonaliśmy mnożenie $a_i$ przez $s$ $10\_000$ razy.

Ze względu na mechanizmy optymalizacji wykonywania programów przez procesor (m.in spekulacyjne wykonywanie instrukcji (*speculative execution*), przwidywanie rozgałęzień programu (*branch predicton*)) musieliśmy zmniejszyć przewidywalność naszego kodu. 
W tym celu kolejne wartości $i$ potrzebne w wielomianie $a_i$ wybieramy z tablicy, która została uprzednio wymieszana. Ponadto po każdym wykonaniu mnożenia losujemy liczbę. Dzięki temu uniemożliwiamy procesorowi przewidywanie dalszych kroków programu.

## Analiza danych

Analiza pliku \verb|.csv|, utworzonego podczas uruchamiania programu w C, polegała na kilku głównych etapach. Po pierwsze - sposób analizy. Nasza analiza skupiała sie głównie na sprawdzaniu, kiedy jest wykonywany tzw. *skok*, czyli wykonanie operacji modulo. Na początku tworzymy więc tabelę z wszystkimi możliwymi współczynnikami oraz kiedy będą przechodzić do tej operacji. Tablica ta ma rozmiar $Q \times Q$, a zapisywane są do niej liczba $1$ jeśli skok został wykonany lub $0$ jesli owego nie było.

Kolejnym krokiem była analiza danych z pliku. Aktualnie, dzięki dużym różnicom w wykonaniu się funkcji z oraz bez modulo, działamy wyłącznie na medianie oraz różnicy median z kolejnych liczb. Dzięki temu możemy zaobserwować skoki.

Następnie utworzyliśmy profil dla naszego analizowanego przypadku. Profilowanie polega na przypisaniu kolejnych wartości różnic do prostej, do której wybrany punkt ma bliżej. Występują dwie proste, odpowiednio $y=min(y)$ oraz $y=max(y)$, a reszta punktów *znajduje się w ich środku*. Dzięki temu otrzymamy profil danego współczynnika.

Ostatnim krokiem jest sprawdzenie, do którego współczynnika z tablicy $Q \times Q$ wyliczony profil ma najbliżej. Dzięki temu otrzymujemy finalnie poprawny tajny współczynnik.
