import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import math
matplotlib.rcParams['figure.figsize'] = (15, 10)

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


def plot_eigenvector(vector, m):
    # create a heatmap plot of eigenvectors
    matplotlib.rcParams['figure.figsize'] = (15, 15)
    plt.imshow(abs(vector),extent=[0,m/2,-m/4,m/4])
    plt.colorbar(orientation='vertical')
    plt.xlabel("Matrix Element, m")
    plt.ylabel("Charge State, n")
    plt.title("Matrix plot of Eigenvectors")
    plt.savefig("figs/matrix_plot_of_eigenvectors.png")
    plt.show()

def plot_eigenvalues(values,ej=52):
    matplotlib.rcParams['figure.figsize'] = (15, 10)
    plt.plot(np.real(values)[:40] / ej)
    plt.xlabel("Matrix Element, m")
    plt.ylabel("Em/Ej")
    plt.savefig("figs/plot_of_eigenvalues.png")
    plt.show()

s = 80
simple_vals, simple_vectors = tridiagonal_matrix(s)
plot_eigenvector(simple_vectors, m=s)
plot_eigenvalues(simple_vals)


# -------------- analytical ----------------------------
def analytical_energy(k, Ec=0.214, Ej=52, ng=0.0):
    q = -2 * Ej / Ec
    r = k + 1 - (k + 1) % 2 + 2 * ng * (-1) ** k
    mat = 4 * Ec * scipy.special.mathieu_a(int(r), q)
    return mat

# this needs more work. I don't think it's working at all. Probably should do in mathematica instead
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

def plot_analytical_energy(p_values, m_values, ej=52):
    #plots mathematica and python analytical values
    matplotlib.rcParams['figure.figsize'] = (15, 10)
    plt.plot(np.real(p_values) / ej, label='Python')
    plt.plot(m_values/ ej, label='Mathematica')
    plt.xlabel("Matrix Element, m")
    plt.ylabel("Em/Ej")
    plt.legend()

ana = np.array([analytical_energy(i, ng=1E-10) for i in np.arange(0, s)])
mathematica = np.loadtxt("data_ng_nearly0.dat", dtype=float)
plot_analytical_energy(ana,mathematica)
plt.savefig("figs/compare_python_mathematica.png")
plt.show()
#^^ the python part of this graph is really weird?
# it's kind of backwards
# maybe don't use this
# compare the python and mathematica functions using ng=0 exactly, comment that this is a special case

ana0 = np.array([analytical_energy(i, ng=0.0) for i in np.arange(0, s)])
mathematica0 = np.loadtxt("data_ng_0.dat", dtype=float)
plot_analytical_energy(ana0,mathematica0)
plt.savefig("figs/compare_python_mathematica_ng0.png")
plt.show()

# compare mathematica and simple values
ej=52
matplotlib.rcParams['figure.figsize'] = (15, 8)
fig,ax=plt.subplots(1,2)
ax[0].plot(simple_vals[:40], label='Numerical')
plt.xlabel("Matrix Element, m")
plt.ylabel("Em/Ej")
plt.legend()
plt.plot(mathematica[:40]/16.0, label='Analytical', color='orange')
plt.xlabel("Matrix Element, m")
plt.ylabel("Em/Ej")
plt.legend()
plt.savefig("figs/compare_analytical_numerical.png")
plt.show()

# ---------------- koch method ----------------

def k_function(m, ng=0.0):
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

#t = np.loadtxt("data2.dat")
#plt.plot(t)
#plt.show()

def koch_wavefunction(m, ng=0.000001,Ec=0.214, Ej=52):
    # finds the wavefunction according to koch paper
    # for phi = -2pi to 2pi`
    #for phi in :
        #psi = np.exp(1j *ng * phi)/np.sqrt(2.0)
        #mat = scipy.special.mathieu_cem(-2*(ng-k_function(m,ng)),-Ej/2*Ec,phi/2.0)[0]
    r = -2*(ng-k_function(m,ng))
    print(m,k_function(m,ng),int(r))
    q = -Ej/2*Ec
    return np.array([(np.exp(1j *ng * phi)/np.sqrt(2.0))*scipy.special.mathieu_cem(int(r),q,phi/2.0)[0] for phi in np.arange(-2*np.pi, 2*np.pi, 0.1)])


for m in np.arange(41,80,1):
    wave_test = koch_wavefunction(m)
    plt.plot(wave_test, label='%s' %m)
plt.legend()
plt.show()

for m in np.arange(41,80,1):
    wave_test = koch_wavefunction(m)
    prob = np.absolute(wave_test ** 2)
    plt.plot(prob, label='%s' %m)
plt.legend()
plt.show()

"""
m=99
wave_test = koch_wavefunction(m)
plt.plot(wave_test, label='%s' %m)
plt.legend()
plt.show()

prob = np.absolute(koch_wavefunction(41))**2
plt.plot(prob)
plt.show()
"""