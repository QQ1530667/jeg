# Packages import
import numpy as np

# Parameter set
## physical constants
m=9.109*10**-31
hb=1.054*10**-34
C=1.602*10**-19
p=0.2988*4
## controlable variables
E=np.arange(0.0,p,0.000001)
V0=np.tile([0.0],len(E))
V1=np.tile([p],len(E)) #height of the potential well
alpha=0.09014 #effective coefficient of mass
beta=0.0657
a=4.0 #nm
b=10.0
y1=np.zeros(len(E))
y2=np.zeros(len(E))
theta=np.zeros(len(E))
phi=np.zeros(len(E))
for i in range(len(E)):
    theta[i]=np.sqrt(2*alpha*m*C*(V1[i]-E[i]))/hb*a*1E-9
    phi[i]=np.sqrt(2*beta*m*C*(E[i]-V0[i]))/hb*b*1E-9
    y1[i]=abs(np.exp(-theta[i])/2-np.cos(phi[i])/(np.sin(phi[i])+1))
    y2[i]=abs(np.exp(-theta[i])/2-np.cos(phi[i])/(np.sin(phi[i])-1))

minarg1=y1.argsort()[:11]
minarg2=y2.argsort()[:11]
for i in minarg1:
    print(E[i],y1[i])
    
for i in minarg2:
    print(E[i],y2[i])
    
    # y1[i]=theta[i]-np.log(2*np.cos(phi[i])/(np.sin(phi[i])+1))
    # y2[i]=theta[i]+np.log(2*np.cos(phi[i])/(np.sin(phi[i])-1))