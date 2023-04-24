# -*- coding: utf-8 -*-
"""
Created on Tue Apr  18 21:13:27 2023

@author: 33641
"""



# The error is truly uncertain!V0 and p and n are difficult to set!V0 better to be small, p be large!
# Packages import
import time
import numpy as np
import matplotlib.pyplot as plt
start=time.time()

# Parameter set
## physical constants
me=9.109*10**-31
hb=1.054*10**-34
C=1.602*10**-19

## controlable variables
V0=0.29988*4 #height of the potential well
# V0=0.19988 #height of the potential well
alpha=0.09014 #effective coefficient of mass
beta=0.0657
a=4.0 #nm
b=10.0
n=5 #number of the rigion
V=np.append(np.tile([0,V0],int((n-1)/2)),[0])        
m=np.insert([1.0,1.0],1,np.append(np.tile([alpha,beta],int((n-3)/2)),[alpha]))*me

# construct the x boundary list
def xFunc():
    xb=np.tile([a,b],int((n-1)/2))
    x=np.zeros(n-1)
    for i in range(1,n-1):
        x[i]=np.sum(xb[0:i])
    return x

# Recruit calculation of amplitude
A=np.zeros((n,2),dtype=complex)
A[n-1][0]=1
AA=A.copy() # temperary storage
B=np.zeros(n-1,dtype=complex)
p=5000  #number of points gathered
# E=np.setdiff1d(np.logspace(-5,np.log(V0+0.2),p),V)  # a very strange thing is that at the turning point, there is a calculating influence of the lowest state transmisson!!!!
E=np.setdiff1d(np.linspace(0,V0+0.001,p),V)  # a very strange thing is that at the turning point, there is a calculating influence of the lowest state transmisson!!!!
k=np.zeros(n,dtype=complex)
T=np.zeros(len(E))
x=xFunc()
K=np.zeros((len(E),n),dtype=complex) #check k whether is symmetrical in this potential!
for l in range(len(E)):
    for i in range(n):
        k[i]=np.sqrt(complex(2*m[i]*C*(E[l]-V[i])))/hb*10**-9
        K[l][i]=k[i] 
    B[n-2]=(m[n-2]*k[n-1])/(m[n-1]*k[n-2])    
    A[n-2][0]=(1+B[n-2])*np.exp(1j*(k[n-1]-k[n-2])*x[n-2])*A[n-1][0]/2
    A[n-2][1]=(1-B[n-2])*np.exp(1j*(k[n-1]+k[n-2])*x[n-2])*A[n-1][0]/2
    AA[n-2][0]=A[n-2][0]
    AA[n-2][1]=A[n-2][1]
    s=n
    while (s-3>=0):
        B[s-3]=(m[s-3]*k[s-2])/(m[s-2]*k[s-3])
        AA[s-3][0]=((1+B[s-3])*np.exp(1j*(k[s-2]*x[s-3]))*AA[s-2][0]+(1-B[s-3])*np.exp(-1j*(k[s-2]*x[s-3]))*AA[s-2][1])*np.exp(-1j*(k[s-3]*x[s-3]))/2
        AA[s-3][1]=((1-B[s-3])*np.exp(1j*(k[s-2]*x[s-3]))*AA[s-2][0]+(1+B[s-3])*np.exp(-1j*(k[s-2]*x[s-3]))*AA[s-2][1])*np.exp(1j*(k[s-3]*x[s-3]))/2        
        A[s-3][0]=AA[s-3][0]
        A[s-3][1]=AA[s-3][1]
        AA[s-2][0]=AA[s-3][0]
        AA[s-2][1]=AA[s-3][1]
        s=s-1
    T[l]=abs(A[n-1][0])**2/abs(A[0][0])**2
# T=np.real(T)
plt.plot(E,T)
M=np.vstack((T,E))
Ma=M[:,np.argsort(-M[0,:])].transpose() #left line is T,right is E
if Ma[0][0]>1:
    print("Something wrong!!!")
else:
    print("successfull!")

end=time.time()
print("running time:%.2f Seconds"%(end-start))
for i in range(1,len(T)-1):
    if E[i]<V0:
        if T[i]>=T[i-1] and T[i]>=T[i+1]:
            print(E[i],T[i])   
# abs(np.sort(-T)) #problem is T=3??

# print(np.stack(E,T))
# calculate the transmision
# print(AFunc(n))
# Draw