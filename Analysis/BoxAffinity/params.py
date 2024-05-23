import sys,os
import numpy as np
from mpmath import fp
from scipy.special import gamma

#--math constants

euler=fp.euler 

#--color factors

CA=3.0
CF=4.0/3.0
TR=0.5
TF=0.5

#--set_masses

me   = 0.000511
mmu  = 0.105658
mtau = 1.77684
mu   = 0.055
md   = 0.055
ms   = 0.2
mc   = 1.51
mb   = 4.92
mZ   = 91.1876
mW   = 80.398
M    = 0.93891897
MD   = 1.875613/2.    
MN   = 0.93956536
Mpi  = 0.135
Mpip =  0.13957039 
Mpim =  0.13957039 
Mpi0 =  0.1349768 
Mk   = 0.497 
Mkp =  0.493677
Mkm =  0.493677
Mk0 =  0.497648 

me2   = me**2 
mmu2  = mmu**2 
mtau2 = mtau**2
mu2   = mu**2  
md2   = md**2  
ms2   = ms**2  
mc2   = mc**2  
mb2   = mb**2  
mZ2   = mZ**2  
mW2   = mW**2  
M2    = M**2  
Mpi2  = Mpi**2

#--electroweak couplings 

c2w = mW2/mZ2
s2w = 1.0-c2w
s2wMZ = 0.23116
alfa  = 1/137.036
alphaSMZ = 0.118
alfa2  = alfa**2

#--quark charges

eU2 = 4.0/9.0
eD2 = 1.0/9.0
couplings={1:eU2,2:eD2,3:eD2,4:eU2,5:eD2,6:eU2}
  
#--RG params

Q0 =1.0
Q20=1.0
alphaS0=0.2  #--input alphaS(Q20) for forward evolution 
aZ=alphaSMZ/(4*np.pi)

#--TMD evolution params
bTcut_ff=1.0
b0   = 2e-1
bmax = 1.0
g0 = 0.1
C1=1.
C2=1.
C4=1.
C5=0.5






















