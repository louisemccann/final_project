import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special

# ---------- numerical -------------------
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


# -------------- analytical ----------------------------
def analytical_energy(k, Ec=0.214, Ej=52, ng=0):
    q = -2*Ej / Ec
    r = k+1 - (k+1)%2 + 2*ng*(-1)**k
    mat = 4*Ec*scipy.special.mathieu_a(r,q)
    return mat


def analytical_wavefunction(k, Ec=0.214, Ej=52, ng=0):
    # generates a wavefunction between -pi and pi
    Ek = analytical_energy(k, Ec, Ej, ng)
    const1 = int(4*Ek/Ec)
    const2 = -2*Ej/Ec
    print(Ek,const1,const2)
    pi_range = np.arange(-np.pi, np.pi,0.01)
    mathieu_c = np.array([scipy.special.mathieu_cem(const1, const2, theta/2.0)[0] for theta in pi_range])
    mathieu_s = np.array([scipy.special.mathieu_sem(const1, const2, theta/2.0)[0] for theta in pi_range])
    ex = np.exp([(1j *ng *theta)/np.sqrt(2*np.pi) for theta in pi_range])
    return ex*(mathieu_c + (1j*(-1)**(k+1))*mathieu_s)
    #return ex, mathieu_c, mathieu_s


a = analytical_wavefunction(1,1,1)
#print(a)
plt.plot(np.arange(-np.pi, np.pi,0.01), abs(a)**2)
plt.show()