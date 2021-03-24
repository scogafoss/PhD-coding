# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 19:59:45 2020

@author: Oli

"""
import numpy as np
import matplotlib.pyplot as plt
l2=[]
steps=[]
logl2 = []
logsteps = []
# read_file = open("l2_normalised_correctx2.txt","r")
read_file = open("l2_data_gem.txt","r")
for line in read_file:
    steps.append((float(line.split()[0])))
    l2.append((float(line.split()[1])))
#print(l2)
logl2=np.log(l2)
logsteps=np.log(steps)
loc2=np.where(logsteps>-5)[0]
loc1=np.where(logsteps<0)[0]
# print(loc1[len(loc1)-1],loc2[0])
tempx=logsteps[loc1[0]:loc2[len(loc2)-1]]
# print(loc2[len(loc2)-1])
# print(logsteps[0:5])
tempy=logl2[loc1[0]:loc2[len(loc2)-1]]
meanx=np.average(tempx)
meany=np.average(tempy)
store1=0
store2=0
for i in range(len(tempx)):
    store1=store1+((tempx[i]-meanx)*(tempy[i]-meany))
    # print(store1)
    store2=store2+(tempx[i]-meanx)**2
    # print(store2)
# print(store1,store2)#gradient=store1/store2
gradient=store1/store2
plt.xlabel("delta x")
plt.ylabel("L2")
#plt.xticks(np.arange(1,200,25))
#plt.yticks(np.arange(0,200,25))
#fig, ax = plt.subplots()
#ax.set_xticks(np.linspace(0,200, 10))
#ax.set_yticks(np.linspace(0,200,10))
plt.plot(steps,l2)
plt.figure()
plt.xlabel("ln(delta x)")
plt.ylabel("ln(L2)")
# plt.plot(logsteps,logl2,linestyle='none',marker='.', markersize='2')
plt.plot(logsteps,logl2,marker='.', markersize='2')
print(gradient)
plt.title('gradient = '+str(gradient))
# plt.title(' steps')