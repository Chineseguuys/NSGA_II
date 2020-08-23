import numpy as np 
from matplotlib import pyplot as plt 


x = np.arange(0.0,1.0,0.01)
y = 1 - np.sqrt(x)
file = open("log.txt")
val_list = file.readlines()

a = val_list[0].split(' ')[:-1]
b = val_list[1].split(' ')[:-1]

a = np.array(a, dtype = float)
b = np.array(b,dtype = float)
print(a)
print(b)

plt.scatter(a,b, marker='*', label = "NSGAII")
#plt.plot(x,y, color = 'red', label = u"理论")
plt.grid(which = 'both')
plt.legend(loc = 'best')
plt.title("ZDT6")
plt.xlabel("f1")
plt.ylabel("f2")
plt.show()



