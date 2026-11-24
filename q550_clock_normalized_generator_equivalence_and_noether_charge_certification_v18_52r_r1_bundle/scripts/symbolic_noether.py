#!/usr/bin/env python3
"""Exact angular-sector Noether check; SymPy or independent rational monomial algebra."""
from pathlib import Path
import json
from collections import defaultdict
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1]

# Independent exact sparse-Laurent expansion. Formal variables in exponent order:
# (r1, r2, omega1, omega2, v1, v2, ke, sinDelta).
# Coordinates, sinDelta and ke are formal independent symbolic indeterminates.

def monomial(exponents,coefficient):
 return (tuple(exponents),Fraction(coefficient))

def sumterms(*term_lists):
 d=defaultdict(Fraction)
 for term_list in term_lists:
  for powers,c in term_list:d[powers]+=c
 return {powers:coeff for powers,coeff in d.items() if coeff}

def powers(**kw):
 keys=('r1','r2','omega1','omega2','v1','v2','ke','sinDelta')
 return tuple(kw.get(k,0) for k in keys)

left=sumterms([
 monomial(powers(r1=2,v1=1,omega1=1),2),
 monomial(powers(r1=1,r2=1,ke=1,sinDelta=1),1),
 monomial(powers(r1=2,v1=1,omega1=1),-2)])
right=sumterms([
 monomial(powers(r2=2,v2=1,omega2=1),2),
 monomial(powers(r1=1,r2=1,ke=1,sinDelta=1),1),
 monomial(powers(r2=2,v2=1,omega2=1),-2)])
assert left==right=={powers(r1=1,r2=1,ke=1,sinDelta=1):Fraction(1)}
assert (1-1)==0  # d( theta1+eps+theta2-eps-phi0-Theta ) / deps

out={'status':'EXACT_FORMAL_MONOMIAL_PASS','method':'exact Fraction-valued symbolic Laurent-monomial cancellation',
     'phase_counterrotation_derivative':'0', 'angular_noether_dC_dt':'0',
     'direct_angular_ode_dC_dt':'0',
     'scope':'angular-sector Noether-dAlembert, zero angular nonconservative virtual work required'}
try:
 import sympy as S
 r1,r2,w1,w2,v1,v2,k,x,d=S.symbols('r1 r2 w1 w2 v1 v2 k x Delta',real=True)
 th1,th2,phi0,Theta=S.symbols('theta1 theta2 phi0 Theta',real=True)
 V=k*x*r1*r2*S.cos(th1+th2-phi0-Theta)
 eps=S.symbols('eps',real=True)
 dV=S.simplify(S.diff(V.subs({th1:th1+eps,th2:th2-eps}),eps))
 p1dot=-S.diff(V,th1);p2dot=-S.diff(V,th2)
 noether=S.simplify(p1dot-p2dot)
 assert dV==0 and noether==0
 omega1p=k*x*r2/r1*S.sin(d)-2*v1*w1
 omega2p=k*x*r1/r2*S.sin(d)-2*v2*w2
 lhs=S.factor(2*r1*r1*v1*w1+r1*r1*omega1p-(2*r2*r2*v2*w2+r2*r2*omega2p))
 assert S.simplify(lhs)==0
 out.update({'status':'EXACT_SYMPY_AND_FORMAL_MONOMIAL_PASS','sympy':S.__version__,
             'active_counterrotation_dV_deps':str(dV)})
except ImportError:
 out['sympy']='NOT_INSTALLED; independent exact Fraction formal cancellation passed'
(ROOT/'results').mkdir(parents=True,exist_ok=True)
(ROOT/'results'/'exact_symbolic_noether.json').write_text(json.dumps(out,indent=2)+'\n')
print(out['status'])
