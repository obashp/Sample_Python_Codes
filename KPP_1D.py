import numpy as np
import matplotlib.pyplot as plt
import torch
from matplotlib import cm


Pi = 4.0*np.arctan(1.0)

def Set_Bc1(V, N, fw, fe, dx):
    V[0] = V[1] - fw*dx
    V[N-1] = V[N-2] + fe*dx
    
def Set_Bc(V,N,ubl,ubr):
    V[0] = 2.0*ubl-V[1]
    V[N-1] = 2.0*ubr-V[N-2]
    
def Compute_Fe(u,alphax,i,dx):
    Fe = alphax[i]*(u[i+1]-u[i])/dx
    return Fe
    
def Flux_Integral(u,alphax,dx,N,FI):
    Fe = 0.0
    for i in range(N-1):
        Fw = -Fe
        Fe = Compute_Fe(u,alphax,i,dx)
        FI[i] = (Fe + Fw)/dx
        
    return FI


def Source(u,r,beta,N,S):
    S = r*u.mul(1-u**beta)
    return S

def SinSource(x,frq,S):
    S = 0.5*(1-torch.cos(2*Pi*frq*x))
    return S
    
def RK2(U0,ubl,ubr,alphax,r,beta,N,dx,ntime,dt,nstore):
    UU = torch.zeros([N,nstore+1],dtype = torch.float32)
    t = torch.zeros([nstore+1],dtype = torch.float32)
    FI = torch.zeros([N],dtype=torch.float32)
    S = torch.zeros([N],dtype=torch.float32)
    storetime = ntime//nstore
    
    UU[:,0] = U0
    U = U0
    tr = 0.0
    
    Set_Bc(U,N,ubl,ubr)
    istr = 1;
    for itime in range(1,ntime+1):
        FI = Flux_Integral(U,alphax,dx,N,FI)
        S = Source(U,r,beta,N,S)
        Ui = U + 0.5*(FI+S)*dt
        Set_Bc(Ui,N,ubl,ubr)
        FI = Flux_Integral(Ui,alphax,dx,N,FI)
        S = Source(Ui,r,beta,N,S)
        U = U + (FI+S)*dt
        Set_Bc(U,N,ubl,ubr)
        tr = tr+dt
        UU[:,itime] = U
        if(itime % 200 == 0):
            print(itime,t[itime].item(),torch.max(U).item())
        if(itime%storetime == 0):
            UU[:,istr] = U
            t[:,istr] = tr
    
    return UU,t


def Sol(X,N,r):
    k = torch.sqrt(torch.tensor(r/6))
    u = (1+torch.exp(k*X))**(-2)
    return u
   
    
def compute_forwardSolution(U0,ubl,ubr,N,dx,dt,Ntime,alphax,r,beta): 
    Uc,tc = RK2(U0,ubl,ubr,alphax,r,beta,N,dx,Ntime,dt)
    
    Uc[0,:] = ubl; Uc[-1,:] = ubr
    return Uc,tc

def Solution_diff(Ue,Uc,N,Ntime):
    L2_diff = torch.zeros(Ntime+1,dtype=torch.float32)
    L2_diff = 1.0/np.sqrt(N-2) * torch.norm(Ue[1:-1,:]-Uc[1:-1,:],dim=0)
        
    return L2_diff

def Get_params(N,alphax,xmin,xmax):
    dx = (xmax-xmin)/(N-2)
    # dt = 0.25*(0.5*dx**2/torch.max(alphax))
    
    Xc = torch.linspace(xmin-0.5*dx,xmax+0.5*dx,N,dtype=torch.float32)
    
    return dx,Xc

def test_diffusion(ngrids):       
    xmin = 0.0; xmax = 1.0

    r = 0.0
    beta = 0.0
    
    
    for igrid in range(ngrids):
        N = 2**(7+igrid)+2
        Ntime = 5000*2**(2*igrid)
        alphax = torch.ones(N,dtype=torch.float32)
        dt, dx, Xc = Get_params(N, alphax, xmin, xmax)
        #Compute exact solution
        Ue = torch.zeros([N,Ntime+1],dtype=torch.float32)
        for itime in range(0,Ntime+1):
            k1 = -alphax*Pi**2*itime*dt
            k2 = 100*k1
            u1 = torch.exp(k1)*torch.sin(Pi*Xc) + 0.1*torch.exp(k2)*torch.sin(10*Pi*Xc)
            Ue[:,itime] = u1

        #Set Boundary conditions
        ubl = 0.0;ubr = 0.0
        Ue[0,:] = ubl; Ue[-1,:] = ubr
        
        #Compute solution using RK2
        U0 = Ue[:,0]
        Uc,tc = compute_forwardSolution(U0, ubl, ubr, N, dx, dt, Ntime, alphax, r, beta)
        Xc[0] = 0.0; Xc[-1] = 1.0

        #Test solution difference
        L2_diff = Solution_diff(Ue, Uc, N, Ntime)
        plt.plot(tc,L2_diff)





    
# Finalize grid size of 256 for forward model.
# Solve backward model on a grid of 320 and 512

#Extrapolate solution
    

    

    
    

        