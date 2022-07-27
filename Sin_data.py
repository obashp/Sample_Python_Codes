import numpy as np
import torch 
import h5py
import matplotlib.pyplot as plt


Pi = 4.0*np.arctan(1.0)

def Set_Bc(V,N):
    V[0] = -V[1]
    V[N-1] = 2.0-V[N-2]
    
def d2U(u,dx,N,FI):
    for i in range(1,N-1):
        FI[i] = (u[i+1] - 2*u[i] + u[i-1])/dx**2
    
    return FI

def SinSource(x,frq,S):
    S = 0.5*(1-torch.cos(2*Pi*frq*x))
    return S

def Source(u,r,beta,N,S):
    S = r*u.mul(1-u**beta)
    return S

def RK2_Sin(U0,alphax,frq,x,N,dx,ntime,dt,nstore):
    UU = torch.zeros([N,nstore+1],dtype = torch.float32)
    t = torch.zeros([nstore+1],dtype = torch.float32)
    FI = torch.zeros([N],dtype=torch.float32)
    S = torch.zeros([N],dtype=torch.float32)
    storetime = (ntime+1)//nstore
    
    UU[:,0] = U0
    U = U0
    tr = 0.0
    
    Set_Bc(U,N)
    istr = 1;
    for itime in range(1,ntime+1):
        FI = d2U(U,dx,N,FI)
        S = SinSource(x, frq, S)
        Ui = U + 0.5*(FI+S)*dt
        Set_Bc(Ui,N)
        FI = d2U(Ui,dx,N,FI)
        S = SinSource(x, frq, S)
        U = U + (FI+S)*dt
        Set_Bc(U,N)
        tr = tr+dt
        if(itime % 200 == 0):
            print(itime,tr,torch.max(U).item())
        if(itime%storetime == 0):
            UU[:,istr] = U
            t[istr] = tr
            istr = istr+1
    
    return UU,t

def RK2(U0,alphax,r,beta,N,dx,ntime,dt,nstore):
    UU = torch.zeros([N,nstore+1],dtype = torch.float32)
    t = torch.zeros([nstore+1],dtype = torch.float32)
    FI = torch.zeros([N],dtype=torch.float32)
    S = torch.zeros([N],dtype=torch.float32)
    storetime = (ntime+1)//nstore
    
    UU[:,0] = U0
    U = U0
    tr = 0.0
    
    Set_Bc(U,N)
    istr = 1;
    for itime in range(1,ntime+1):
        FI = d2U(U,dx,N,FI)
        S = Source(U,r,beta,N, S)
        Ui = U + 0.5*(FI+S)*dt
        Set_Bc(Ui,N)
        FI = d2U(Ui,dx,N,FI)
        S = Source(Ui,r,beta,N, S)
        U = U + (FI+S)*dt
        Set_Bc(U,N)
        tr = tr+dt
        if(itime % 200 == 0):
            print(itime,tr,torch.max(U).item())
        if(itime%storetime == 0):
            UU[:,istr] = U
            t[istr] = tr
            istr = istr+1
    
    return UU,t

N = 1026
alphax = torch.ones(N,dtype=torch.float32)
xmin = 0; xmax = 1
ubl = 0.0;ubr = 0.0
dx = (xmax - xmin)/(N-2)
Xc = torch.linspace(xmin-0.5*dx,xmax+0.5*dx,N,dtype=torch.float32)
Xc[0] = xmin; Xc[-1] = xmax

U0 = torch.zeros(N,dtype=torch.float32)
r = 100
c1 = np.sqrt(r/6)
wv = 5*np.sqrt(r/6)
for i in range(N):
    if(Xc[i] > 0.5-1.0/16 and Xc[i] < 0.5 + 1.0/16):
        U0[i] = torch.sin(4*Pi*(Xc[i]+1.0/16))**2
    elif(Xc[i] >= 0.5 + 1.0/16):
        U0[i] = 1.0


dt = 0.4*0.5*dx**2
frq = 1
tf = 1/16
ntime = int(tf//dt)+1
nstore = 8

beta = 1


Uf, t = RK2(U0, alphax, r, beta, N, dx, ntime, dt, nstore)

Ue = torch.zeros([N,nstore+1],dtype=torch.float32)
Ue[:,0] = U0
for i in range(1,nstore+1):
    X = Xc - wv*t[i]
    Ue[:,i] = (1 + torch.exp(c1*X))**(-2)

fname = 'Solution_KPP_C2_%d.h5'%N
hf = h5py.File(fname, 'w')
for i in range(nstore+1):
    datastr = 'dataset_%d' % i
    hf.create_dataset(datastr, data=Uf[:,i])
hf.create_dataset('time',data=t)
hf.create_dataset('Xc',data=Xc)
hf.close()