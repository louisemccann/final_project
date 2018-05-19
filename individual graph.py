import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import scipy.sparse
import math

s=80
def analytical_energy(k, Ec=0.214, Ej=52, ng=0.0):
    q = -2 * Ej / Ec
    r = k + 1 - (k + 1) % 2 + 2 * ng * (-1) ** k
    mat = 4 * Ec * scipy.special.mathieu_a(int(r), q)
    return mat


def plot_analytical_energy(p_values, m_values, ej=52):
    # plots mathematica and python analytical values
    #matplotlib.rcParams['figure.figsize'] = (15, 10)
    plt.plot(np.real(p_values) / ej, label='Python')
    plt.plot(m_values / ej, label='Mathematica')
    plt.xlabel("Energy State, m")
    plt.ylabel("$E_m / E_J$")
    plt.legend()


ana = np.array([analytical_energy(i, ng=0) for i in np.arange(0, s)])
mathematica = np.loadtxt("data_ng_0.dat", dtype=float)
plot_analytical_energy(ana, mathematica)
plt.savefig("figs/comparepythonmathematica.png")
plt.show()