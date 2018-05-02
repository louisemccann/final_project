import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.special

# load data from mathematica

ej = 52.0
ec=52.0
mathematica = np.loadtxt("data2.dat", dtype=float)
matplotlib.rcParams['figure.figsize'] = (15, 10)
plt.plot(mathematica / ej, label='Mathematica')
plt.savefig("figs/mathematica_koch.png")
plt.show()

def koch_wavefunction(m, ng=0.000001, Ec=0.214, Ej=52):
    # finds the wavefunction according to koch paper
    # for phi = -2pi to 2pi
    # for phi in :
    # psi = np.exp(1j *ng * phi)/np.sqrt(2.0)
    # mat = scipy.special.mathieu_cem(-2*(ng-k_function(m,ng)),-Ej/2*Ec,phi/2.0)[0]
    r = -2 * (ng - k_function(m, ng))
    #print(m, k_function(m, ng), int(r))
    q = -Ej / 2 * Ec
    return np.array(
        [(np.exp(1j * ng * phi) / np.sqrt(2.0)) * scipy.special.mathieu_cem(int(abs(r)), q, phi / 2.0)[0] for phi in
         np.arange(-2 * np.pi, 2 * np.pi, 0.1)])


ng=0
for m in range(2):
    r = -2 * (ng - mathematica[m])
    q = -ej / 2 * ec
    print(r, q)
    wave = np.array(
        [(np.exp(1j * ng * phi) / np.sqrt(2.0)) * scipy.special.mathieu_cem(int(abs(r)), q, phi / 2.0)[0] for phi in
         np.arange(-2 * np.pi, 2 * np.pi, 0.1)])
    plt.plot(wave)
    plt.show()


mathematica = np.loadtxt("data3.dat", dtype=float)
matplotlib.rcParams['figure.figsize'] = (15, 10)
plt.plot(mathematica / ej, label='Mathematica')
plt.savefig("figs/mathematica_koch.png")
plt.show()
