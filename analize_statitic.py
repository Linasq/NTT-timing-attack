import pandas as pd
import os
import numpy as np


def find_jumps_for_given_Q(Q, refresh = 0):
    if not os.path.isfile(f"./Q_tables/{Q}_Q.csv") or refresh == 1:
        with open(f"./Q_tables/{Q}_Q.csv", "wb") as f:
            for a in range(Q):
                a_copy = 2*a
                for j in range(Q-1):
                    if a_copy > Q:
                        f.write(b'1,')
                        a_copy -= Q
                    else:
                        f.write(b'0,')
                    a_copy += a
                f.seek(-1, os.SEEK_END)
                f.write(b'\n')

    jumps = pd.read_csv(f"./Q_tables/{Q}_Q.csv", header=None)
    return jumps

def get_differences_from_c_data(filename):
    df = pd.read_csv('statistic_cycle_measurement.csv', header=None)
    cols = df.columns
    for col in cols:
        df[col] = df[col].astype(float)

    average = df.median(axis=1)
    print(average)
    differences = np.diff(average.to_numpy())
    return differences, 0, max(differences), average

def create_jumps_profile(averages, minimum, maximum):
    profile = []
    for a in averages:
        if maximum - a < a - minimum:
            profile.append(1)
        else:
            profile.append(0)
    return profile

pd.options.display.max_rows = 999
pd.options.display.max_columns = 999
Q = 97
jumps_df = find_jumps_for_given_Q(Q)
differences, min_diff, max_diff, av = get_differences_from_c_data('statistic_cycle_measurement.csv')
print(differences, min_diff, max_diff)
profile = create_jumps_profile(differences, min_diff, max_diff)
print(profile)
result = jumps_df[jumps_df.eq(profile).all(axis=1)]
print(result)

import matplotlib.pyplot as plt

plt.scatter(range(len(av)), av, s=90)

plt.xlabel("Numer pomiaru")
plt.ylabel("Średnia wartość")
plt.title("Wykres średnich wartości")
plt.grid(True)
plt.show()


plt.scatter(range(len(differences)), differences, s=90)

plt.xlabel("Numer pomiaru")
plt.ylabel("Średnia różnica")
plt.title("Różnice między kolejnymi pomiarami")
plt.grid(True)
plt.show()