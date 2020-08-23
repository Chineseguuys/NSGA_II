import numpy as np 
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


ax = plt.subplot(111, projection='3d')
file = open("log.txt")
val_list = file.readlines()

a = val_list[0].split(' ')[:-1]
b = val_list[1].split(' ')[:-1]
c = val_list[2].split(' ')[:-1]

a = np.array(a,dtype = float)
b = np.array(b, dtype = float)
c = np.array(c, dtype = float)

ax.scatter(a,b,c,marker = '.')
plt.title("DTLZ2")
plt.show()
