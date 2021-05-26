import matplotlib.pyplot as plt
# import numpy as np
read=open('phi_2d.txt','r')
phi=[]
x=[]
for line in read:
    phi.append(float(line.split()[0]))
    x.append(float(len(phi)-1))
plt.plot(x,phi)