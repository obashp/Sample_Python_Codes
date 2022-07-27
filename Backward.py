import KPP_1D as Kp
import numpy as np
import torch 
import torch.nn as nn
from torch.autograd import grad
import torch.nn.functional as F
import matplotlib.pyplot as plt
import scipy
from scipy import sparse
from scipy.sparse import coo_matrix
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import h5py

def interp_Solution(U0,Xinit,Xfinal):
    f = interpolate.interp1d(Xinit,U0,kind='cubic')
    Ufinal = f(Xfinal)
    return torch.tensor(Ufinal,dtype=torch.float32)

def Add_Noise(U0,N,percent_noise):
    WN = torch.tensor(np.random.normal(0,percent_noise/100,N))
    U0 = U0+WN
    return U0

def Solve_TDMA(A,b,N):
    #forward substitution
    for i in range(N-1):
        A[i][2] = A[i][2]/A[i][1]
        b[i] = b[i]/A[i][1]
        A[i][1] = 1.0
        
        A[i+1][1] = A[i+1][1] - A[i][2]*A[i+1][0]
        b[i+1] = b[i+1] - b[i]*A[i+1][0]
        A[i+1][0] = 0.0
        
    b[N-1] = b[N-1]/A[N-1][1]
    A[N-1][1] = 1.0
    #backward substitution
    for i in range(N-2,-1,-1):
        b[i] = b[i] - A[i+1][2]*b[i+1]
        A[i+1][2] = 0.0
    
    return b

class Backward_Solution:
    def __init__(self,Nb,xmin,xmax,ubl,ubr,alphax,initpt,finalpt,fname,isource=0):
        self.xmin = xmin
        self.xmax = xmax
        self.N = Nb
        self.dx = (xmax - xmin)/(Nb-2)
        self.source = isource
        self.Ntime = 8*(Nb-2)
        
        self.X = torch.linspace(xmin-0.5*self.dx,xmax+0.5*self.dx,Nb)
        self.X[0] = xmin; self.X[-1] = xmax
        self.Uinv = torch.zeros(Nb,dtype=torch.float32)
        self.FI = torch.zeros(Nb,dtype=torch.float32)
        self.S = torch.tensor([])
        if(self.source != 0):
            self.S = torch.zeros(self.N,dtype=torch.float32)
            
        self.Init_Solution(initpt, finalpt, fname)
        self.Init_Bc(ubl, ubr)
        self.alpha = alphax
            
        
    def Init_Bc(self,ubl,ubr):
        self.ubl = ubl
        self.ubr = ubr
        
    def Set_Bc(self,U):
        U[0] = 2*self.ubl - U[1]
        U[-1] = 2*self.ubr - U[-2]
        
    def Set_SinSourcefreq(self,freq):
        if(self.source == 1):
            self.frq = freq
        else:
            self.frq = 0.0
        
    def Set_SourceParams(self,r=1,beta=1):
        self.r = r
        self.beta = beta        
        
    def Set_Rp(self,Rp):
        self.Rp = Rp
        
    def Init_Solution(self,initpt,finalpt,fname):
        hf = h5py.File(fname,'r')
        initstr = 'dataset_%d' % initpt
        U0 = hf.get(initstr)
        U0 = torch.tensor(U0)
        initstr = 'dataset_%d' % finalpt
        Uf = hf.get(initstr)
        Uf = torch.tensor(Uf)
        Xc = hf.get('Xc')
        Xc = torch.tensor(Xc)
        tf = hf.get('time')
        tf = torch.tensor(tf)
        t0 = tf[initpt]
        t1 = tf[finalpt]
        hf.close()
        self.dt = torch.abs(t0-t1)/self.Ntime
        self.U0inv = interp_Solution(U0,Xc,self.X)
        self.Ufinv = interp_Solution(Uf,Xc,self.X)
    
    def Add_Noise(self,noise):
        self.U0inv = Add_Noise(self.U0inv,self.N,noise)
        
    def Compute_FI(self,U):
        self.FI = Kp.Flux_Integral(U, self.alpha, self.dx, self.N, self.FI)
        
    def Compute_Source(self,U):
        if(self.source == 1):
            self.S = Kp.SinSource(self.X,self.frq,self.S)
        elif(self.source == 2):
            self.S = Kp.Source(U, self.r, self.beta, self.N, self.S)
            
    def Time_Advance(self,U,Ntime):
        A = torch.zeros([self.N,3],dtype=torch.float32)
        b = torch.zeros(self.N, dtype=torch.float32)
        
        r = self.alpha*self.dt/self.dx**2
        for itime in range(1,Ntime+1):
            A[:,0] = -0.5*r
            A[:,1] = 1 + 0.5*2.0*r
            A[:,2] = -0.5*r
            A[0,0] = 0.0; A[0,1] = 1.0; A[0,2] = 1.0
            A[-1,2] = 0.0; A[-1,1] = 1.0; A[-1,0] = 1.0
            
            self.Compute_FI(U)
            self.Compute_Source(U)
            b = self.FI
            if(self.source != 0):
                b = b + self.S
            dU = Solve_TDMA(A, b*self.dt, self.N)
            U = U+dU
            self.Set_Bc(U)
            if(itime%200 == 0):
                print(itime,torch.max(U).item())
        
        return U
    
    def Solve_Backward(self):
        KU = self.U0inv.clone()
        Uinv = self.U0inv.clone()
        KU = self.Time_Advance(KU, self.Ntime)
        for iCG in range(10):
            print("K2U")
            K2U = self.Time_Advance(Uinv, 2*self.Ntime)
            K2U = K2U + self.Rp*Uinv
            r0 = KU - K2U
            if iCG == 0:
                nr0 = torch.dot(r0,r0).item()
            s0 = r0.clone()
            
            #Compute a0
            print("K2s")
            K2s = self.Time_Advance(s0, 2*self.Ntime)
            K2s = K2s + self.Rp*s0
            a0 = torch.dot(r0,r0)/torch.dot(K2s,s0)
            print(a0.item())
            #Compute Next iteratie
            Uinv = Uinv + a0*s0
            r1 = r0 - a0*K2s
            
            b0 = torch.dot(r1,r1)/torch.dot(r0,r0)
            s1 = r1 + b0*s0
            
            nr1 = torch.dot(r1,r1).item()
            loss = torch.sqrt(F.mse_loss(a0*s0,0*s0)/len(s0)).item()
            print(iCG,nr1,nr1/nr0,loss)
            r0 = r1; s0 = s1;
            if(nr1/nr0 < 1e-7):
                break
        
        self.Uinv = Uinv
    
    def Error(self,errfilename):
        den = F.mse_loss(self.Uinv,0*self.Uinv)
        misfit = self.Time_Advance(self.Uinv, self.Ntime)
        num = F.mse_loss(misfit,self.U0inv)
        
        #Save the data
        with open(errfilename,'a') as f:
            fstr = '%12.10lf\t%12.10lf\t%12.10lf\n'%(self.Rp,num,den)
            f.write(fstr)
    
        error = self.Uinv - self.Ufinv
        return error
    
    def Save_solution(self,filename,errfname):
        error = self.Error(errfname)
        hf = h5py.File(filename, 'w')
        hf.create_dataset('Beta',data=self.Rp)
        hf.create_dataset('Uinv',data=self.Uinv)
        hf.create_dataset('error',data=error)
        hf.create_dataset('X',data=self.X)
        hf.close()
    

    
Nb = 130
xmin = 0.0; xmax = 1.0
ubl =  0.0; ubr =  0.0
alphax = torch.ones(Nb,dtype=torch.float32)
initpt = 4; finalpt = 2
fname = 'Solution_Sin_4.h5'
frq = 4
noise = 2.0

BS = Backward_Solution(Nb, xmin, xmax, ubl, ubr, alphax, initpt, finalpt, fname,1)
BS.Set_SinSourcefreq(frq)
# BS.Set_SourceParams()
BS.Add_Noise(noise)
BS.Set_Rp(BS.dx**2)
BS.Solve_Backward()

fname = 'Inv_data_Sin_noise_%d_%dx%d.h5' % (frq,initpt,finalpt)
errfname = 'misfit_solnorm.dat'
BS.Save_solution(fname, errfname)


