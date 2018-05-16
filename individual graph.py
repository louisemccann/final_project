import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import scipy.sparse
import math
matplotlib.rcParams['figure.figsize'] = (15, 10)
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
s = 200
simple_vals, simple_vectors = tridiagonal_matrix(s)
f=[]
for m in range(200-1):
    n = np.matrix(simple_vectors[m]).T
    charge = np.diag(range(-100,100))
    n1 = np.matrix(simple_vectors[m+1])
    f1=np.matmul(np.matmul(n1,charge),n)
    f.append(abs(f1.flat[0]))
plt.plot(f[100:110], label='Numerical Result')
plt.ylabel("Coupling Constant")
plt.xlabel("Energy State, m")
#plt.savefig("figs/coupling.png")
#plt.show()

def koch(s):
    return [0.5*(math.sqrt(j+1/2.0)*(52/8*0.214)**0.25)+2.7 for j in np.arange(0,s,0.1)]

plt.plot(np.arange(0,10,0.1),koch(10), label='Approximation')
plt.legend()
plt.savefig('figs/couplingcomparekoch.png', Transparent=True)
plt.show()