import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special


def tridiagonal_matrix(x, Ec=0.214, Ej=52, ng=0):
    # where x is the size of the matrix
    # x should be even
    # allows for negative values of n
    diag = [4 * Ec * (m - ng) ** 2 for m in range(int(-x / 2), int(x / 2))]
    off_diag = [-Ej / 2] * (x - 1)
    vals, vectors = scipy.linalg.eigh_tridiagonal(diag, off_diag)
    vals = np.real(vals)
    return vals, vectors

def plot_eigenvector(vector):
    # create a heatmap plot of eigenvectors
    plt.imshow(abs(vector))
    plt.colorbar(orientation='vertical')
    plt.xlabel("m")
    plt.ylabel("n")
    plt.title("Matrix plot of Eigenvectors")
    plt.show()


s = 100
matplotlib.rcParams['figure.figsize'] = (15, 8)
simple_vals, simple_vectors = tridiagonal_matrix(s)
plot_eigenvector(simple_vectors)
