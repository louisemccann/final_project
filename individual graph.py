import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import scipy.sparse
import math

matplotlib.rcParams['figure.figsize'] = (20, 8)
matplotlib.rcParams.update({'font.size': 22})


# ---------- numerical -------------------
def tridiagonal_matrix(x, Ec=0.214, Ej=52, ng=0.0):
    # where x is the size of the matrix
    # x should be even
    # allows for negative values of n
    diag = [4 * Ec * (m - ng) ** 2 for m in range(int(-x / 2), int(x / 2))]
    off_diag = [-Ej / 2] * (x - 1)
    vals, vectors = scipy.linalg.eigh_tridiagonal(diag, off_diag)
    vals = np.real(vals)
    return vals, vectors

def plot_eigenvalues(values, ej=52):
    #matplotlib.rcParams['figure.figsize'] = (15, 10)
    plt.plot(np.real(values)[:40] / ej)
    plt.xlabel("Matrix Element, m")
    plt.ylabel("$E_m / E_J$")
    plt.savefig("figs/plot_of_eigenvalues.png")
    plt.show()


s = 80
simple_vals, simple_vectors = tridiagonal_matrix(s)
mathematica0 = np.loadtxt("data_ng_nearly0.dat", dtype=float)

# compare mathematica and simple values
ej = 52
fig, ax = plt.subplots(1, 2)
ax[0].plot(simple_vals[:40]/ej, label='Numerical')
ax[0].set_xlabel("Matrix Element, m")
ax[0].set_ylabel("$E_m / E_J$")
ax[0].legend()
ax[1].plot(mathematica0[:40] / ej, label='Analytical', color='orange')
ax[1].set_xlabel("Matrix Element, m")
ax[1].set_ylabel("$E_m / E_J$")
ax[1].legend()
plt.savefig("figs/compare_analytical_numerical.png")
plt.show()
