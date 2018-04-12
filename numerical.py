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
