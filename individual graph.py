import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import scipy.sparse
import math

root = "C:/Users/Louise/Documents/project/final_project/data/data_ng_"
add = ['0','02']
for i in add:
    inp = root+i+'.dat'
    l = np.loadtxt(inp, dtype=float)
    plt.plot(l/52.0, label='$n_g$ = 0.' +i)
plt.legend()
plt.show()