import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib
import scipy.special
import scipy.sparse
import math

matplotlib.rcParams['figure.figsize'] = (15, 15)
matplotlib.rcParams.update({'font.size': 22})

root = 'data/data_ng_0'
add = ['',1,2,3,4,5]
for i in add:
    plt.plot(np.loadtxt(root+str(i)+'.dat')/52, label=str(i))
plt.legend()
plt.show()

