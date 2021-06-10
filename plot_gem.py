# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 21:18:04 2020

@author: Oli
"""


import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
x=[]
phi=[]
phi2=[]
xgem=[]
phigem = []
phigem2=[]
factor_test=[]
difference_test=[]

read_file = open("phi_2d.txt","r")
for line in read_file:
    x.append((float(line.split()[2])))
    phi.append((float(line.split()[0])))
    phi2.append((float(line.split()[1])))
read_file.close()
read=open("xgem_ygem1_ygem2.txt","r")
for line in read:
    xgem.append((float(line.split()[0])))
    phigem.append((float(line.split()[1])))
    phigem2.append((float(line.split()[2])))
read.close()
#print(l2)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.xlabel("x")
plt.ylabel("phi")
# plt.xticks(np.arange(1,200,25))
# plt.yticks(np.arange(0,200,25))
#fig, ax = plt.subplots()
#ax.set_xticks(np.linspace(0,200, 10))
#ax.set_yticks(np.linspace(0,200,10))
plt.plot(x,phi,marker='x',label='Code1')
plt.plot(xgem,phigem,label='Gem1',marker='o')
plt.plot(x,phi2,label='Code2',marker='x')
plt.plot(xgem,phigem2,label='Gem2',marker='o')
# import matplotlib.lines as mlines
# print(phi)
# blue_line = mlines.Line2D(xgem,phigem, color='blue',markersize=15, label='Gem')
# red_line = mlines.Line2D(x,phi, color='red',markersize=15, label='FE')
# orange_line = mlines.Line2D(xgem,phigem2, color='orange',markersize=15, label='Gem2')
# green_line = mlines.Line2D(x,phi2, color='green',markersize=15, label='FE2')
# plt.legend(handles=[blue_line,red_line,orange_line,green_line])
# plt.legend(handles=[blue_line,red_line])
plt.legend()
plt.title('Code vs Gem,nodes='+str(len(x)))
plt.show()
factor_test.append(phi[0]/phigem[0])
factor_test.append(phi[len(phi)-1]/phigem[len(phigem)-1])
print('factors= ',factor_test)
difference_test.append(phi[0]-phigem[0])
difference_test.append(phi[len(phi)-1]-phigem[len(phigem)-1])
print('differences= ',difference_test)