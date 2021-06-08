import matplotlib.pyplot as plt
# import numpy as np
read=open('phi_2d.txt','r')
phi=[]
x=[]
xgem=[0,0.2,0.4,0.6,.8,1,1.2,1.4,1.6,1.8,2.0]
#vacuum
# ygem=[ 2.2117159E-01, 3.2021168E-01,3.8039529E-01, 4.2921357E-01,4.8812101E-01,5.8652427E-01, 6.7500091E-01,7.0000478E-01, 6.7672206E-01,5.9551096E-01,4.1802028E-01]
#reflective
# ygem = [1.0856936E+00,1.0961847E+00,1.1302268E+00,1.1961549E+00,1.3101119E+00,1.5000000E+00,1.6898881E+00,1.8038451E+00,1.8697732E+00,1.9038153E+00,1.9143064E+00]
#tb reflective lr vacuum
ygem =[2.7585347E-01, 3.9894183E-01,4.7154538E-01,5.2993410E-01, 6.0327674E-01,7.2821248E-01,8.4226402E-01,8.7751668E-01,8.5158136E-01,7.5150176E-01, 5.2728193E-01]
for line in read:
    phi.append(float(line.split()[0]))
    x.append(float(line.split()[1]))
plt.plot(x,phi,label='code')
plt.plot(xgem,ygem,label='gem')
plt.legend()

# import matplotlib.pyplot as plt
# # import numpy as np
# read=open('phi_2d.txt','r')
# phi=[]
# x=[]
# xgem=[]
# ygem=[]
# for line in read:
#     phi.append(float(line.split()[0]))
#     x.append(float(line.split()[1]))
#     xgem.append(float(line.split()[2]))
#     ygem.append(float(line.split()[3]))
# plt.plot(x,phi,label='code')
# plt.plot(xgem,ygem,label='gem')
# plt.legend()