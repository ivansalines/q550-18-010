#!/usr/bin/env python3
"""Fast algebraic tests; no Q550 physical ODE integration."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import metric_path_attribution as audit
import numpy as np
class TestGeometry(unittest.TestCase):
    def test_three_factor_shapley_exact(self):
        rng=np.random.default_rng(1550)
        basis=rng.normal(size=(3,45));interaction=rng.normal(size=45)
        vertices={}
        for c in (0,1):
            for w in (0,1):
                for m in (0,1):
                    key=frozenset(k for k,x in [('center',c),('width',w),('metric',m)] if x)
                    vertices[key]=basis[0]*c+basis[1]*w+basis[2]*m+interaction*c*w*m*.2
        d,tot,cl=audit.four_interactions_shapley(vertices,('center','width','metric'))
        self.assertLess(cl,1e-13)
        self.assertLess(np.linalg.norm(d['center']+d['width']+d['metric']-tot),1e-13)
        for i,name in enumerate(('center','width','metric')):
            self.assertLess(np.linalg.norm(d[name]-basis[i]-interaction*.2/3),1e-13)
    def test_two_factor_shapley_exact(self):
        rng=np.random.default_rng(44);a,b,c=(rng.normal(size=10) for _ in range(3))
        v={frozenset():np.zeros(10),frozenset(('center',)):a,
           frozenset(('metric',)):b,frozenset(('center','metric')):a+b+c}
        d,total,cl=audit.four_interactions_shapley(v,('center','metric'))
        self.assertLess(cl,1e-13)
        self.assertLess(np.linalg.norm(d['center']-(a+c/2)),1e-13)
    def test_kinetic_exterior_covariance(self):
        rng=np.random.default_rng(28);q,_=np.linalg.qr(rng.normal(size=(5,5)))
        B=rng.normal(size=10);A=rng.normal(size=(5,5));G=A.T@A+np.eye(5)
        L=audit.mg.wedge_map(q)
        H=audit.mg.induced(G)
        G1=q@G@q.T;H1=audit.mg.induced(G1)
        self.assertAlmostEqual(float(B@H@B),float((L@B)@H1@(L@B)),places=9)
    def test_angles_not_additive(self):
        a=np.array([1.,0.]);b=np.array([0.,1.]);c=np.array([-1.,0.])
        self.assertAlmostEqual(audit.scalar_angle(a,b),90)
        self.assertAlmostEqual(audit.scalar_angle(b,c),90)
        # The signed vector decomposition, not angle addition, carries the identity.
    def test_oriented_identity(self):
        X=np.array([[1.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                    [.9,.1,0,0,0,0,0,0,0,0],
                    [.8,.2,0,0,0,0,0,0,0,0]])
        arc,P,u=audit.mg.oriented_polygon(X,np.eye(10))
        self.assertTrue(np.isfinite(arc));self.assertEqual(P.shape,(45,))
if __name__=='__main__':unittest.main()
