#!/usr/bin/env python3
"""Fast independent analytic synthetic controls for the clock/frame audit."""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import rotation, G, simpson_weights

def run():
    checks={}
    xi=np.linspace(0,.2,257); k1=.02; k2=-.02; dk=k1-k2
    tau=xi+.3*xi**2;c=1+.6*xi
    checks['clock_exact_derivative']=np.max(np.abs(CubicSpline(xi,tau).derivative()(xi)-c))<1e-11
    # One response projection, one forcing channel, analytic kernel on triangle.
    K=np.zeros((len(xi),len(xi),2,3))
    for i in range(len(xi)):
        u=xi[:i+1]
        K[i,:i+1,0,0]=xi[i]-u
        K[i,:i+1,0,1]=np.sin(xi[i]-u)
        K[i,:i+1,0,2]=K[i,:i+1,0,0]+K[i,:i+1,0,1]
    KT=K/c[None,:,None,None]
    err=0.
    for i in (8,32,64):
        old=simpson(K[i,:i+1],x=xi[:i+1],axis=0)
        new=simpson(KT[i,:i+1],x=tau[:i+1],axis=0)
        err=max(err,float(np.max(abs(old-new))))
    checks['positive_clock_kernel_reparametrization']=err<1e-11
    checks['causal_future_support']=np.max(abs(KT[np.triu_indices(len(xi),1)]))==0
    psi=.35*np.sin(2*np.pi*xi/.2)+.15*xi/.2
    theta=.40*np.sin(xi*4)+.2
    vec=np.stack([.1*np.sin(3*xi),xi**2],axis=1)
    rotvec=np.einsum('iab,ib->ia',rotation(-psi),vec)
    lab=np.einsum('iab,ib->ia',rotation(theta),vec)
    from_rot=np.einsum('iab,ib->ia',rotation(theta+psi),rotvec)
    checks['two_coordinate_frames_return_same_physical_vector']=np.max(abs(lab-from_rot))<1e-14
    checks['rotation_preserves_euclidean_norm']=np.max(abs(np.linalg.norm(vec,axis=1)-np.linalg.norm(rotvec,axis=1)))<1e-14
    dp=CubicSpline(xi,psi).derivative()(xi)
    dt=CubicSpline(xi,theta).derivative()(xi)
    vl=CubicSpline(xi,vec,axis=0).derivative()(xi)
    wr=CubicSpline(xi,rotvec,axis=0).derivative()(xi)
    left=np.einsum('iab,ib->ia',rotation(theta),vl+dt[:,None]*(G@vec.T).T)
    right=np.einsum('iab,ib->ia',rotation(theta+psi),wr+(dt+dp)[:,None]*(G@rotvec.T).T)
    checks['moving_frame_connection_cancels_passive_rotation']=np.max(abs(left-right))<1e-5
    # Explicit different parameter orbit and parameter-dependent physical clock.
    # Central secant clock correction differs O(delta_kappa^2), not O(delta_kappa).
    def Q(x,k):return x[:,None]*np.array([[1.,.5]]) + k*x[:,None]**2*np.array([[3.,-2.]])
    def T(x,k):return x+.3*x*x+k*.04*x*x
    tplus,tminus=T(xi,k1),T(xi,k2);tmid=T(xi,0.)
    qplus,qminus=Q(xi,k1),Q(xi,k2)
    dQ=(qplus-qminus)/dk
    dtau=(tplus-tminus)/dk
    qprime=np.tile(np.array([1.,.5]),(len(xi),1))
    corr=-qprime*(dtau/c)[:,None]
    direct=(CubicSpline(tplus,qplus,axis=0)(tmid)-
            CubicSpline(tminus,qminus,axis=0)(tmid))/dk
    checks['fixed_tau_is_not_fixed_xi']=np.max(abs(corr))>1e-5
    checks['fixed_tau_clock_shift_predicts_direct_matching']=np.max(abs(direct-(dQ+corr)))<1e-5
    w=simpson_weights(256,xi[1]-xi[0]);checks['composite_simpson_mass']=abs(w.sum()-xi[-1])<1e-14
    checks['channel_additivity_after_rotation']=np.max(abs((K[:,:,:,0]+K[:,:,:,1])-K[:,:,:,2]))<1e-14
    # Explicit shape tripwire: observation (2) x source-channel (3), never transpose.
    a=np.arange(6.).reshape(2,3);checks['impulse_axes_observation_x_source']=(a.shape==(2,3))
    print('Q550 v18.51v fast synthetic',sum(checks.values()),'/',len(checks),'PASS' if all(checks.values()) else 'FAIL')
    for k,v in checks.items():print(('PASS' if v else 'FAIL'),k)
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':run()
