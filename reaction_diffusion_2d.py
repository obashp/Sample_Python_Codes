# -*- coding: utf-8 -*-
"""Reaction_Diffusion_2D.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ats8XLE_Qv9scBeVgpgcUvdi_t6niyUm
"""

import io
import numpy as np
import matplotlib.pyplot as plt
import torch
from matplotlib import cm
from torch.autograd import grad
import pandas as pd
from google.colab import files 
import h5py

Pi = 4.0*np.arctan(1.0)

def compute_L2norm(V,Nx,Ny):
    L2 = 0.0
    for icol in range(1,Nx-1):
        for jrow in range(1,Ny-1):
            L2 = L2 + V[icol][jrow]**2
            
    L2 = np.sqrt(L2/((Nx-2)*(Ny-2)))
    return L2

def Set_Bc(V, Nx, Ny):
    for icol in range(1,Nx-1):
        V[icol][0] = -V[icol][1]
        V[icol][-1] = -V[icol][-2]
        
    for jrow in range(1,Ny-1):
        V[0][jrow] = -V[1][jrow]
        V[-1][jrow] = -V[-2][jrow]

def diffusion_coeff_x(x):
    return 1-x**2
    # return 1e-5

def diffusion_coeff_y(y):
    # return 1-y**2
    return 0.01*(1.0-y**2)

def Compute_Fe(u,alphax,icol,jrow,dx):
    Fe = alphax[icol]*(u[icol+1][jrow]-u[icol][jrow])/dx
    return Fe

def Compute_Gn(u,alphay,icol,jrow,dy):
    Gn = alphay[jrow]*(u[icol][jrow+1]-u[icol][jrow])/dy
    return Gn

def Flux_Integral(u,alphax, alphay, dx,dy,Nx,Ny,FI):
    Fe = np.zeros(Ny,dtype=float)
    icol = 0; jrow = 1;
    while(icol < Nx-1):
        jrow = 1
        while(jrow < Ny-1):
            Fw = -Fe[jrow]
            Fe[jrow] = Compute_Fe(u,alphax,icol,jrow,dx)
            FI[icol][jrow] = (Fe[jrow]+Fw)/dx
            jrow = jrow+1
        icol = icol+1
    
    Gn = np.zeros(Nx,dtype=float)
    icol = 1; jrow = 0;
    while(jrow < Ny-1):
        icol = 1
        while(icol < Nx-1):
            Gs = -Gn[icol]
            Gn[icol] = Compute_Gn(u,alphay,icol,jrow,dy)
            FI[icol][jrow] = FI[icol][jrow] + (Gn[icol]+Gs)/dy
            icol = icol+1
        jrow = jrow+1
        
    return FI

def Source(u,r,beta,Nx,Ny,S):
  S = r*u.mul(1-u**beta)
    # for icol in range(1,Nx-1):
    #     for jrow in range(1,Ny-1):
    #         S[icol][jrow] = r*u[icol][jrow]*(1-u[icol][jrow]**beta)
  return S

def RK2(U0,alphax,alphay,r,beta,Nx,Ny,dx,dy,ntime,dt):
    N = (Nx)*(Ny)
    UU = torch.zeros([N,ntime+1],dtype = torch.float32)
    t = torch.zeros([ntime+1],dtype = torch.float32)
    FI = torch.zeros([Nx,Ny],dtype=torch.float32)
    S = torch.zeros([Nx,Ny],dtype=torch.float32)
    
    
    UU[:,0] = U0.reshape(N)
    U = U0
    Set_Bc(U,Nx,Ny)
    
    for itime in range(1,ntime+1):
        
        FI = Flux_Integral(U,alphax,alphay, dx, dy, Nx, Ny, FI)
        S = Source(U,r,beta,Nx,Ny,S)
        Ui = U + 0.5*(FI+S)*dt
        Set_Bc(Ui,Nx,Ny)
        FI = Flux_Integral(Ui,alphax,alphay,dx, dy, Nx, Ny, FI)
        S = Source(Ui,r,beta,Nx,Ny,S)
        U = U + (FI+S)*dt
        Set_Bc(U,Nx,Ny)
        UU[:,itime] = U.reshape(N)
        if(itime % 50 == 0):
          print(itime, t[itime-1].item()+dt, torch.max(U).item(),torch.min(U).item())
        t[itime] = t[itime-1]+dt
        
    return UU,t

def Add_Noise(u,N,percent_noise):
    WN = torch.tensor(np.random.normal(0,percent_noise/100,N))
    u = u+WN
    return u

Nx = 52; Ny = 52
dx = 1.0/(Nx-2); dy = 1.0/(Ny-2)
alphax = 0.1*torch.ones(Nx,dtype=torch.float32)
alphay = 0.1*torch.ones(Ny,dtype=torch.float32)
Xc = torch.zeros(Nx,dtype=torch.float32)
Yc = torch.zeros(Ny,dtype=torch.float32)

for icol in range(Nx):
    Xc[icol] = (icol-0.5)*dx
    # alphax[icol] = diffusion_coeff_x(icol*dx)
    
for jrow in range(Ny):
    Yc[jrow] = (jrow-0.5)*dy
    # alphay[jrow] = diffusion_coeff_y(jrow*dy)


u = torch.zeros([Nx,Ny],dtype=torch.float32)
for icol in range(Nx):
    for jrow in range(Ny):
        x = Xc[icol]
        y = Yc[jrow]
        # u[icol][jrow] = 1.0*torch.exp(-10*((x-0.5)**2 + (y-0.5)**2))
        u[icol][jrow] = torch.sin(Pi*x)*torch.abs(torch.sin(2*Pi*y))


X,Y = np.meshgrid(Xc.numpy(),Yc.numpy())
plt.contour(X,Y,u)
# print(torch.max())
dt = 2e-4

Uc, tc = RK2(u,alphax,alphay,1.5,2.0,Nx,Ny,dx,dy,2000,dt)

u = Uc[:,800].reshape(Nx,Ny)
plt.contour(X,Y,u,levels=[0,0.01,0.1,0.2,0.4,0.5])
torch.max(u)

hf = h5py.File('Uc1.h5', 'w')
hf.create_dataset('gridsize',data=[Nx,Ny])
hf.create_dataset('dataset_1', data=Uc.detach().numpy())
hf.close()

