# Harmonic oscillator shooting method solution for eigenvalues & eigenfunctions
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.integrate import simps
from scipy.optimize import brentq
L=10.0
h=0.01
x=np.arange(-L,L,h)
e=np.arange(0.95,10.1,0.0001)
all_zeros=[]
# Function for determining eigenenergy of harmonic oscillator in shooting method
def psiAtL(e):
    po,p1=0.0,1.0
    x=-L+h
    while x<=+L:
        po,p1=p1,(2-h**2*(e-x**2))*p1-po
        x+=h
    return p1
# plotting given function
plt.plot(e,psiAtL(e))
plt.grid(True)
plt.show()
# First five eigenenergies evaluated
# we need to find the zeros for the symmetrical case
s=np.sign(psiAtL(e))
for i in range(len(s)-1):
    if s[i]+s[i+1]==0:
        zero=brentq(psiAtL,e[i],e[i+1])
        all_zeros.append(zero)
print(all_zeros)
# RK4 for solving wave function

psi=0;psid=1. # initial psi and first derivative of psi
I=np.array([psi,psid])
xi=-4.;xf=+4.;h1=0.001;x1=xi
e=all_zeros[0]
xs=[]; psi_s=[]
def Rk4(x1,I):
    p=dI(x1,I)
    q=dI(x1+h1/2,I+h1/2*p)
    r=dI(x1+h1/2,I+h1/2*q)
    s=dI(x1+h1,I+h1*r)
    I+=(p+2*q+2*r+s)/6*h1
    x1+=h1
    return x1,I
def dI(x1,I):
    I=psi,psid
    dpsi_dx=psid
    dpsid_dx=-(e-x1**2)*psi
    return np.array([dpsi_dx,dpsid_dx])
while x1<xf:
    x1,I=Rk4(x1,I)
    psi,psid=I
    xs.append(x1)
    psi_s.append(psi)
arr_xs=np.asarray(xs)
arr_psi=np.asarray(psi_s)
sq_arr_psi=arr_psi**2
# normalization constant using simpsons rule from scipy
print("scipy based simpsons Integration Result as normalization constant =",simps(sq_arr_psi,dx=h1))
# Normalization constant
N=1./simps(sq_arr_psi,dx=h1)
arr_psi=N*arr_psi
plt.plot(arr_xs,arr_psi)
plt.show()
    

























    




























