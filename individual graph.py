import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import scipy.sparse
import math

matplotlib.rcParams['figure.figsize'] = (15, 10)
matplotlib.rcParams.update({'font.size': 30})


def tridiagonal_matrix(x, Ec=0.214, Ej=52, ng=0.0):
    # where x is the size of the matrix
    # x should be even
    # allows for negative values of n
    diag = [4 * Ec * (m - ng) ** 2 for m in range(int(-x / 2), int(x / 2))]
    off_diag = [-Ej / 2] * (x - 1)
    vals, vectors = scipy.linalg.eigh_tridiagonal(diag, off_diag)
    vals = np.real(vals)
    return vals, vectors


def plot_eigenvector(vector, m, sub,fig):
    # create a heatmap plot of eigenvectors
    matplotlib.rcParams['figure.figsize'] = (15, 15)
    im=sub.imshow(abs(vector), extent=[0, m / 2, -m / 2, m / 2],aspect='auto')
    fig.colorbar(im,orientation='vertical')
    sub.xlabel("Matrix Element, m")
    sub.ylabel("Charge State, n")
    #plt.savefig("figs/matrixplotofeigenvectors.png",transparent=True)
    #plt.show()

def plot_eigenvalues(values,sub, ej=52):
    matplotlib.rcParams['figure.figsize'] = (15, 10)
    sub.plot(np.real(values)[:40] / ej)
    sub.xlabel("Matrix Element, m")
    sub.ylabel("$E_m / E_J$")
    #plt.savefig("figs/plotofeigenvalues.png",transparent=True)
    #plt.show()

fig, (ax1,ax2) = plt.subplots(1,2)
s = 80
simple_vals, simple_vectors = tridiagonal_matrix(s)
plot_eigenvector(simple_vectors, m=s, sub=ax1,fig=fig)
plot_eigenvalues(simple_vals,sub=ax2)
plt.savefig("figs/numerical_test.png",transparent=True)
plt.show()
