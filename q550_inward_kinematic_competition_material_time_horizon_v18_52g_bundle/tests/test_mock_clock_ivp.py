#!/usr/bin/env python3
"""Fast 10D mock linear oscillator verifying genuine clock/IVP/event pathway."""
from pathlib import Path
import importlib.util,sys
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
spec=importlib.util.spec_from_file_location('g52mock',p)
g=importlib.util.module_from_spec(spec);sys.modules[spec.name]=g;spec.loader.exec_module(g)

class DummyJet:
    def structural_jet(self,z,kappa,params,slope):
        y=np.asarray(z[:10]);out=np.zeros(10)
        out[g.QPOS]=y[g.UVEL]
        out[7]=-y[6]  # harmonic oscillator in Cartesian component 3
        return out,None,None

class DummyEngine:
    jet=DummyJet()
    params=type('Params',(),{'m_x':1.})()
    slope=None

class DummyKin:
    def embedding_and_jac(self,y,mx,tau_ref):
        h=np.r_[y[g.QPOS],y[g.UVEL]]
        j=np.zeros((10,10))
        for i in range(5):j[i,g.QPOS[i]]=1;j[i+5,g.UVEL[i]]=1
        return h,j

ctx={'event':0.,'theta':lambda t:0.}
y0=np.zeros(10);y0[1]=1.
y1=y0.copy();y1[7]=.01
eng=DummyEngine();kin=DummyKin()
a=g.integrate(ctx,eng,0.,kin,1.,y0,2.5,2e-11,2e-12,.02)
b=g.integrate(ctx,eng,0.,kin,1.,y1,2.5,2e-11,2e-12,.02)
times=np.linspace(0,2.5,611)
eval=lambda t:g.compare_at_tau(t,a,b,ctx,eng,0.,kin,1.)
vals=[eval(t) for t in times]
ev=g.find_events(times,vals,eval)
assert len(ev['convexity_zero'])==2,ev
assert abs(ev['convexity_zero'][0]['tau']-np.pi/4)<1e-8
assert abs(ev['convexity_zero'][1]['tau']-3*np.pi/4)<1e-8
assert len(ev['position_separation_maxima'])==1,ev
assert abs(ev['position_separation_maxima'][0]['tau']-np.pi/2)<1e-8
assert abs(vals[0]['D'])<1e-25
assert max(z['material_identity_position_error'] for z in vals)<1e-12
assert abs(vals[50]['delta_xi'])<1e-12
print('PASS: mock authentic-style 10D+material-clock IVP, dense roots C/G/Dmax')
