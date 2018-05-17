import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import scipy.sparse
import math

matplotlib.rcParams['figure.figsize'] = (15, 10)
matplotlib.rcParams.update({'font.size': 22})


def tridiagonal_matrix(x, Ec=0.214, Ej=52, ng=0.0):
    # where x is the size of the matrix
    # x should be even
    # allows for negative values of n
    diag = [4 * Ec * (m - ng) ** 2 for m in range(int(-x / 2), int(x / 2))]
    off_diag = [-Ej / 2] * (x - 1)
    vals, vectors = scipy.linalg.eigh_tridiagonal(diag, off_diag)
    vals = np.real(vals)
    return vals, vectors


def plot_eigenvector(vector, m):
    # create a heatmap plot of eigenvectors
    # matplotlib.rcParams['figure.figsize'] = (15, 15)
    plt.imshow(abs(vector), extent=[0, m / 2, -m / 2, m / 2], aspect='auto')
    plt.colorbar(orientation='vertical')
    plt.xlabel("Matrix Element, m")
    plt.ylabel("Charge State, n")
    plt.savefig("figs/matrixplotofeigenvectors.png", Transparent=True)
    plt.show()


def plot_eigenvalues(values, ej=52):
    # matplotlib.rcParams['figure.figsize'] = (15, 10)
    plt.plot(np.real(values)[:40] / ej)
    plt.xlabel("Matrix Element, m")
    plt.ylabel("$E_m / E_J$")
    plt.savefig("figs/plot_of_eigenvalues.png")
    plt.show()


s = 80
simple_vals, simple_vectors = tridiagonal_matrix(s)
plot_eigenvector(simple_vectors, m=s)
