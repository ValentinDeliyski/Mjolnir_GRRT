import numpy as np
from scipy.optimize import root_scalar

import matplotlib.pyplot as plt

r_ISCO = np.zeros((5000,5000))
M_index = 0
fig1 = plt.figure()
fig1 = fig1.add_subplot(111)
# fig2 = plt.figure()
# fig2 = fig2.add_subplot(111)

COLORS = ["b", "k"]

for M in [1e2, 1e4]:

    compactness_list = np.linspace(1e-4, 1e-0, 5000)

    i = 0

    for compactness in compactness_list:

        a_0 = M / compactness

        def f(r):

            ksi = 2 * a_0 - M + 4
            Y = np.sqrt(M / ksi) * (2 * np.arctan((r + a_0 + M) / np.sqrt(M * ksi)) - np.pi)

            return (1 - 2 / r) * np.exp(Y)

        def dr_f(r):

            ksi = 2 * a_0 - M + 4
            Y = np.sqrt(M / ksi) * (2 * np.arctan((r + a_0 + M) / np.sqrt(M * ksi)) - np.pi)
            dr_Y = 2 / ksi / (1 + (r + a_0 + M)**2 / M / ksi) 

            return 2 / r**2 * np.exp(Y) + f(r) * dr_Y

        def d2r_f(r):

            ksi = 2 * a_0 - M + 4
            Y = np.sqrt(M / ksi) * (2 * np.arctan((r + a_0 + M) / np.sqrt(M * ksi)) - np.pi)
            dr_Y = 2 / ksi / (1 + (r + a_0 + M)**2 / M / ksi) 
            d2r_Y = - 2 / ksi / (1 + (r + a_0 + M)**2 / M / ksi)**2 * 2 * (r + a_0 + M) / M / ksi

            return 2 / r**2 * np.exp(Y) * dr_Y - 4 / r**3 * np.exp(Y) + dr_f(r) * dr_Y + f(r) * d2r_Y
        
        def Isco_condition(r):

            return d2r_f(r) + 3 * dr_f(r) / r - 2 * dr_f(r)**2 / f(r)

        r_ISCO[M_index][i] = root_scalar(Isco_condition, x0 = 6, x1 = 5).root

        i = i + 1

    fit = np.polyfit(compactness_list, r_ISCO[M_index], 10)
    poly = np.poly1d(fit)

    fit_plot = []

    for x in compactness_list:
        fit_plot.append(poly(x))

    fig1.plot(compactness_list, fit_plot, COLORS[M_index])
    fig1.plot(compactness_list, r_ISCO[M_index],'r')
    print(fit)

    M_index = M_index + 1

print(poly(0.5))

plt.show()