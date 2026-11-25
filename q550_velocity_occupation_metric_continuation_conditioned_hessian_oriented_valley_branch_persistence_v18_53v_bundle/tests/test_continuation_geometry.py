#!/usr/bin/env python3
"""Algebraic and synthetic regression tests; no historical IVP execution."""
from pathlib import Path
import numpy as np,sys
from scipy.interpolate import CubicHermiteSpline
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from metric_continuation import analytic_metric,ResidualObjective,mg,improve_root

def checkclose(x,y,tol=3e-7):
    if not np.allclose(x,y,atol=tol,rtol=tol):
        raise AssertionError(f'SYNTHETIC_ERROR {np.max(np.abs(np.asarray(x)-np.asarray(y)))}')
def main():
    S=np.diag([.1,.25,.3,.15,.2]);lam=1.6
    G,Gp,H,Hp=analytic_metric(S,lam)
    h=1e-5
    Gplus,_,Hplus,_=analytic_metric(S,lam+h)
    Gminus,_,Hminus,_=analytic_metric(S,lam-h)
    checkclose(Gp,(Gplus-Gminus)/(2*h),tol=3e-9)
    checkclose(H,(mg.induced(G)),tol=1e-13)
    checkclose(Hp,(Hplus-Hminus)/(2*h),tol=4e-9)
    # Polynomial cumulative source with exact spline representation; the
    # polynomial is not assumed to produce a physical Q550 return.
    t=np.linspace(0,8,1001)
    def cc(x):
        x=np.asarray(x)
        return np.stack([(.08*(i+1))*x +(.007*(i%3+1))*x**2 + .0004*(i+1)*x**3 for i in range(10)],axis=-1)
    def gg(x):
        x=np.asarray(x)
        return np.stack([.08*(i+1)+.014*(i%3+1)*x+.0012*(i+1)*x**2 for i in range(10)],axis=-1)
    prox=CubicHermiteSpline(t,cc(t),gg(t),axis=0)
    obj=ResidualObjective(prox,H,Hp)
    p=np.array([1.5,.45]);f,grad,Hess=obj.value_grad_hess(p)
    a,b=1e-5,1e-4
    numeric_grad=np.array([(obj.value_grad_hess(p+np.eye(2)[j]*a)[0]-obj.value_grad_hess(p-np.eye(2)[j]*a)[0])/(2*a) for j in range(2)])
    numeric_hess=np.column_stack([(obj.value_grad_hess(p+np.eye(2)[j]*b)[1]-obj.value_grad_hess(p-np.eye(2)[j]*b)[1])/(2*b) for j in range(2)])
    checkclose(grad,numeric_grad,2e-8)
    checkclose(Hess,numeric_hess,4e-6)
    fl,gl=obj.lambda_derivatives(p)
    op=ResidualObjective(prox,analytic_metric(S,lam+h)[2],np.zeros((10,10)))
    om=ResidualObjective(prox,analytic_metric(S,lam-h)[2],np.zeros((10,10)))
    checkclose(fl,(op.value_grad_hess(p)[0]-om.value_grad_hess(p)[0])/(2*h),3e-8)
    checkclose(gl,(op.value_grad_hess(p)[1]-om.value_grad_hess(p)[1])/(2*h),5e-8)
    assert np.linalg.eigvalsh(H)[0]>0
    assert abs(np.linalg.det(G)-1)<2e-12
    # Synthetic well-conditioned target minimum: independent quadratic
    # affine residual detects sign error in the implicit response.
    A=np.array([[1.,0.2],[0.3,2.],[0.1,.7]])
    base=np.array([.3,-.5,.2]);pert=np.array([1.,.2,-.3])
    def target(lmb):
        W=np.diag([1.,1.+lmb/10,1.+lmb/4])
        return -np.linalg.solve(A.T@W@A,A.T@W@(base+lmb*pert))
    q=target(1.5)
    xp=(target(1.500001)-target(1.499999))/.000002
    W=np.diag([1.,1+.15,1+.375]);Wp=np.diag([0,.1,.25])
    predicted=-np.linalg.solve(A.T@W@A,A.T@Wp@(A@q+base+1.5*pert)+A.T@W@pert)
    checkclose(predicted,xp,1e-8)
    print('SYNTHETIC_METRIC_DERIVATIVE_HESSIAN_CONTINUATION_PASS')
if __name__=='__main__':main()
