import numpy as np
from scipy.interpolate import CubicSpline

PAIRS5=tuple((i,j) for i in range(5) for j in range(i+1,5))
PAIRS10=tuple((i,j) for i in range(10) for j in range(i+1,10))

def wedge5(A):
    A=np.asarray(A,float)
    return np.array([[A[i,k]*A[j,l]-A[i,l]*A[j,k] for k,l in PAIRS5] for i,j in PAIRS5])

def exterior_generator(A):
    A=np.asarray(A,float)
    L=np.zeros((10,10))
    for n,(k,l) in enumerate(PAIRS5):
        E=np.zeros((5,5));E[k,l]=1.;E[l,k]=-1.
        M=A@E+E@A.T
        L[:,n]=[M[i,j] for i,j in PAIRS5]
    return L

def minrot(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    a=a/np.linalg.norm(a);b=b/np.linalg.norm(b)
    dot=float(a@b)
    if dot<=-0.99:raise RuntimeError('ANTIPODAL_TANGENT_MINIMAL_ROTATION_UNDEFINED')
    K=np.outer(b,a)-np.outer(a,b)
    return np.eye(len(a))+K+(K@K)/(1+dot)

class KinematicFrame:
    def __init__(self,t,v,stride=1):
        t=np.asarray(t,float)[::stride];v=np.asarray(v,float)[::stride]
        self.t=t;self.v=v;self.spline=CubicSpline(t,v,axis=0,extrapolate=False)
        speed=np.linalg.norm(v,axis=1)
        if np.min(speed)<=1e-6:raise RuntimeError('NEAR_ZERO_5D_VELOCITY')
        e=v/speed[:,None]
        self.e=e
        R=np.empty((len(t),5,5));R[0]=np.eye(5)
        # accumulated shortest-turn transport, no imposed torsion in the moving tangent
        for i in range(1,len(t)):
            T=minrot(e[i-1],e[i]);R[i]=T@R[i-1]
        self.R=R
        self.orthogonality=float(np.max(np.abs(np.einsum('nji,njk->nik',R,R)-np.eye(5))))
        self.tangent_error=float(np.max(np.linalg.norm(np.einsum('nij,j->ni',R,e[0])-e,axis=1)))
        if self.orthogonality>1e-8 or self.tangent_error>1e-8:raise RuntimeError('FRAME_DRIFT')
    def at(self,t):
        t=np.asarray(t,float)
        sh=t.shape;tt=t.ravel()
        if tt.min()<self.t[0]-1e-10 or tt.max()>self.t[-1]+1e-10:raise RuntimeError('FRAME_TIME_OUTSIDE_PHYSICAL_TRACE')
        inds=np.clip(np.searchsorted(self.t,tt,side='right')-1,0,len(self.t)-1)
        vv=self.spline(tt)
        ee=vv/np.linalg.norm(vv,axis=1)[:,None]
        ans=np.empty((len(tt),5,5))
        for i in range(len(tt)):ans[i]=minrot(self.e[inds[i]],ee[i])@self.R[inds[i]]
        return ans.reshape(sh+(5,5))
    def omega(self,t):
        v=self.spline(t);vd=self.spline.derivative()(t)
        speed=np.linalg.norm(v)
        e=v/speed; ed=(vd-e*(e@vd))/speed
        return np.outer(ed,e)-np.outer(e,ed)

def rot5(i,j,a):
    c,s=np.cos(a),np.sin(a)
    R=np.eye(5);R[i,i]=R[j,j]=c;R[i,j]=-s;R[j,i]=s
    return R

def rot5_d(i,j,a):
    c,s=np.cos(a),np.sin(a)
    R=np.zeros((5,5));R[i,i]=R[j,j]=-s;R[i,j]=-c;R[j,i]=c
    return R

def gauge(t,t0,t1,A=.22):
    s=(t-t0)/(t1-t0);k=2*np.pi/(t1-t0)
    a=A*np.sin(2*np.pi*s);b=.5*A*np.sin(4*np.pi*s)
    da=A*k*np.cos(2*np.pi*s);db=A*k*np.cos(4*np.pi*s)
    R1=rot5(0,4,a);R2=rot5(1,3,b)
    Rp=rot5_d(0,4,a)*da@R2+R1@rot5_d(1,3,b)*db
    return R1@R2,Rp

def geodesic_polygon(N):
    u=np.asarray(N,float);nn=np.linalg.norm(u,axis=1)
    if nn.min()<=1e-12:raise RuntimeError('ZERO_RESIDUAL_DIRECTION')
    u=u/nn[:,None]
    a,b=u[:-1],u[1:]
    ang=2*np.arctan2(np.linalg.norm(b-a,axis=1),np.linalg.norm(a+b,axis=1))
    k=np.where(ang<1e-5,1+ang**2/6+7*ang**4/360,ang/np.maximum(np.sin(ang),1e-50))
    poly=np.array([np.sum(k*(a[:,i]*b[:,j]-a[:,j]*b[:,i])) for i,j in PAIRS10])
    return float(ang.sum()),poly,u

def angle_deg(a,b):
    aa=np.linalg.norm(a);bb=np.linalg.norm(b)
    if aa<1e-14 or bb<1e-14:return np.nan
    return float(np.degrees(np.arccos(np.clip(a@b/(aa*bb),-1,1))))

if __name__=='__main__':
    rng=np.random.default_rng(11)
    for _ in range(30):
        a=rng.normal(size=5);a/=np.linalg.norm(a);b=a+.1*rng.normal(size=5);b/=np.linalg.norm(b)
        R=minrot(a,b);assert np.linalg.norm(R@a-b)<2e-14
        assert np.linalg.norm(R.T@R-np.eye(5))<2e-14
        A=rng.normal(size=(5,5));A=A-A.T
        g=exterior_generator(A)
        eps=1e-7
        from scipy.linalg import expm
        assert np.linalg.norm((wedge5(expm(eps*A))-np.eye(10))/eps-g)<1e-5
        assert np.linalg.norm(g.T+g)<1e-14
    print('MATH SYNTHETIC PASS')
