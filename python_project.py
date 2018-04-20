import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import math


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


def plot_eigenvector(vector):
    # create a heatmap plot of eigenvectors
    plt.imshow(abs(vector))
    plt.colorbar(orientation='vertical')
    plt.xlabel("m")
    plt.ylabel("n")
    plt.title("Matrix plot of Eigenvectors")
    plt.show()


s = 80
matplotlib.rcParams['figure.figsize'] = (15, 8)
_, simple_vectors = tridiagonal_matrix(s)


#plot_eigenvector(simple_vectors)


# -------------- analytical ----------------------------
def analytical_energy(k, Ec=0.214, Ej=52, ng=0.0):
    q = -2 * Ej / Ec
    r = k + 1 - (k + 1) % 2 + 2 * ng * (-1) ** k
    mat = 4 * Ec * scipy.special.mathieu_a(int(r), q)
    return mat


def analytical_wavefunction(k, Ec=0.214, Ej=52, ng=0.0):
    # generates a wavefunction between -pi and pi
    Ek = analytical_energy(k, Ec, Ej, ng)
    const1 = int(4 * Ek / Ec)
    const2 = -2 * Ej / Ec
    pi_range = np.arange(-np.pi, np.pi, 0.01)
    mathieu_c = np.array([scipy.special.mathieu_cem(const1, const2, theta / 2.0)[0] for theta in pi_range])
    mathieu_s = np.array([scipy.special.mathieu_sem(const1, const2, theta / 2.0)[0] for theta in pi_range])
    ex = np.exp([(1j * ng * theta) / np.sqrt(2 * np.pi) for theta in pi_range])
    return ex * (mathieu_c + (1j * (-1) ** (k + 1)) * mathieu_s)
    # return ex, mathieu_c, mathieu_s


a = analytical_wavefunction(1, 1, 1)
# print(a)
# plt.plot(np.arange(-np.pi, np.pi, 0.01), abs(a) ** 2)
# plt.show()

# ------------------------- compare methods ---------------------
ej = 52.0
simple_vals, _ = tridiagonal_matrix(s, ng=0.0)
simple_vals= simple_vals/25.0
plt.plot(np.real(simple_vals) / ej,label='simple vals')
ana = np.array([analytical_energy(i, ng=0.0) for i in range(0, s)])/100.0
plt.plot(ana/ej, label='analytical')


mathematica = np.loadtxt("data.dat", dtype=float)


plt.plot(mathematica/ej, label='mathematica')
#plt.savefig("test2.png",format='png')
plt.legend()
plt.show()

# ---------------- koch method ----------------

def k_function(m, ng=0.0):
    k = 0
    l = -1.0
    k = round((2 * ng + l / 2.0) % 2) * (round(ng + (l * (-1) ** m) * ((m + 1) // 2.0)))
    l = 1.0
    k += round((2 * ng + l / 2.0) % 2) * (int(ng + (l * (-1) ** m) * ((m + 1) // 2.0)))
    return k


def em_koch(m, ng=0.0,Ec=0.214, Ej=52):
    return Ec * scipy.special.mathieu_a(2, ng + np.real(k_function(m,ng))) * (-Ej/(2*Ec))


k_test = [k_function(m) for m in range(1,40)]
koch_test = [em_koch(m) for m in range(1,40)]
plt.plot(range(1,40),np.array(koch_test)/52.0, label='koch')
#plt.savefig("test.png",format='png')
#plt.show()

t = np.loadtxt("data2.dat")
plt.plot(t)
plt.show()

def koch_wavefunction(m, ng=0.0,Ec=0.214, Ej=52):
    # finds the wavefunction according to koch paper
    # for phi = -2pi to 2pi
    return [(np.exp(1j *ng * phi)/np.sqrt(2.0)) * scipy.special.mathieu_cem(-2*(ng-k_function(m,ng)),-Ej/2*Ec,phi/2.0) for phi in np.arange(-2*np.pi, 2*np.pi, 0.1)]

wave_test = koch_wavefunction(1)
