import matplotlib.pyplot as plt 
import numpy as np

for i in range(750):
    filename='t'+str(i)
    x,w,V=np.loadtxt(filename+'.dat', delimiter=' ', unpack=True)
    plt.xlabel('X')
    plt.plot(x,w, label='Probability Density')
    plt.plot(x,V, label='Potential (Not to Scale)')
    plt.legend()
    plt.savefig(str(i)+'.png')
    plt.close()

