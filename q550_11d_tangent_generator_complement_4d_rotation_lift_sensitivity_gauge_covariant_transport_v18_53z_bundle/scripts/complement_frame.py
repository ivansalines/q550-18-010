#!/usr/bin/env python3
"""SO(5) transport induced by a model-conditioned complementary skew field.

R_y is the authenticated minimal-turn frame from v18.53y. With a choice
of 11D horizontal lift, R_z=R_y U and U'=R_y.T B R_y U.
B annihilates the observed kinetic velocity versor. The lift and tangent
horizon are diagnostic choices; no unique material connection is claimed.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import expm
from connection_math import wedge5

def skew(A):return .5*(A-A.T)

def cayley(A):
    I=np.eye(A.shape[-1]);return np.linalg.solve(I-.5*A,I+.5*A)

class ComplementFrame:
    def __init__(self,t,B,base):
        self.t=np.asarray(t,float); self.base=base
        self.spline=CubicSpline(self.t,np.asarray(B,float),axis=0,extrapolate=False)
        self.U=np.empty((len(t),5,5));self.U[0]=np.eye(5)
        for i in range(len(t)-1):
            mid=.5*(t[i]+t[i+1]);dt=t[i+1]-t[i]
            R=base.at(mid); b=skew(self.spline(mid))
            K=skew(R.T@b@R)
            self.U[i+1]=expm(dt*K)@self.U[i]
        self.max_orthogonality=float(np.max(np.abs(np.einsum('nji,njk->nik',self.U,self.U)-np.eye(5))))
        if self.max_orthogonality>1e-8:raise RuntimeError('COMPLEMENT_FRAME_NOT_ORTHOGONAL')
    def at(self,times):
        times=np.asarray(times,float);flat=times.ravel()
        if flat.min()<self.t[0]-1e-10 or flat.max()>self.t[-1]+1e-10:
            raise RuntimeError('COMPLEMENT_FRAME_OUT_OF_TRACE')
        ind=np.searchsorted(self.t,flat,side='right')-1
        ind=np.clip(ind,0,len(self.t)-2);dt=flat-self.t[ind]
        mid=(self.t[ind]+flat)/2
        Rm=self.base.at(mid)
        Bm=self.spline(mid);Bm=.5*(Bm-np.swapaxes(Bm,-1,-2))
        A=np.einsum('nji,njk,nkl->nil',Rm,Bm,Rm)
        A=.5*(A-np.swapaxes(A,-1,-2))
        M=cayley(dt[:,None,None]*A)
        U=np.einsum('nij,njk->nik',M,self.U[ind])
        return np.einsum('nij,njk->nik',self.base.at(flat),U).reshape(times.shape+(5,5))

def pullback(N,t,frame):
    R=frame.at(t)
    return np.array([wedge5(a.T)@n for a,n in zip(R,N)])
