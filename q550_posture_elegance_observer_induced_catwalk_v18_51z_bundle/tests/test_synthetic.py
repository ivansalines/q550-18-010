#!/usr/bin/env python3
"""Independent synthetic identities for Q550 v18.51z, no historical data."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from run_audit import rot,cross,angles,observer_phase

def verify():
    n=1001;x=np.linspace(0,1,n)
    # Two clockwise channels, increasing omega2 amplitude. Their sum can turn CCW.
    aa=-.4*x;bb=np.pi/2-.4*x
    a=np.stack([np.cos(aa),np.sin(aa)],axis=1)
    b=(.1+1.8*x)[:,None]*np.stack([np.cos(bb),np.sin(bb)],axis=1)
    o=np.zeros_like(a);v=np.stack([a,b,o,a+b],axis=2)
    ang=angles(v)
    tests={}
    tests['both_synthetic_channels_turn_clockwise']=bool(ang[-1,0]<ang[1,0] and ang[-1,1]<ang[1,1])
    tests['synthetic_total_turns_counterclockwise']=bool(ang[-1,3]>ang[1,3])
    tests['synthetic_three_channel_closure']=bool(np.max(abs(v[:,:,0]+v[:,:,1]+v[:,:,2]-v[:,:,3]))<1e-14)
    phi=.31*x+.24*np.sin(4*np.pi*x)
    w=rot(v,phi);wa=angles(w)
    tests['rotated_lab_reconstruction']=bool(np.max(abs(rot(w,-phi)-v))<1e-14)
    tests['oriented_wedge_covariant']=bool(np.max(abs(cross(v[:,:,0],v[:,:,1])-cross(w[:,:,0],w[:,:,1])))<1e-14)
    tests['mutual_angle_covariant']=bool(np.max(abs((ang[1:,1]-ang[1:,0])-(wa[1:,1]-wa[1:,0])))<1e-12)
    da=np.gradient(ang[1:],x[1:],axis=0,edge_order=2)
    dw=np.gradient(wa[1:],x[1:],axis=0,edge_order=2)
    dphi=np.gradient(phi[1:],x[1:],edge_order=2)
    good=[0,1,3]
    tests['rotating_connection_restores_angular_rate']=bool(np.max(abs(dw[:,good]+dphi[:,None]-da[:,good]))<1e-10)
    tests['omitted_connection_creates_defect']=bool(np.max(abs(dw[:,good]-da[:,good]))>.1)
    dv=np.gradient(v,x,axis=0,edge_order=2)
    own=np.stack([cross(v[:,:,j],dv[:,:,j]) for j in range(3)],axis=1)
    pair=sum((cross(v[:,:,i],dv[:,:,j])+cross(v[:,:,j],dv[:,:,i]) for i,j in ((0,1),(0,2),(1,2))))
    total=cross(v[:,:,3],dv[:,:,3]);
    tests['wedge_turning_budget']=bool(np.max(abs(total-own.sum(axis=1)-pair))<1e-12)
    tests['positive_clock_reparam_preserves_sign']=bool(np.all(np.sign(np.gradient(ang[1:,3],x[1:]))==np.sign(np.gradient(ang[1:,3],2*x[1:]))))
    # Reflection is not a proper SO(2) rotation: it reverses oriented area.
    r=v.copy();r[:,1,:]*=-1
    tests['reflection_reverses_oriented_area']=bool(np.max(abs(cross(v[:,:,0],v[:,:,1])+cross(r[:,:,0],r[:,:,1])))<1e-13)
    # Angles at the exact zero vector must not be fabricated.
    zero=v.copy();zero[0,:,:]=0
    tests['zero_response_angle_masked']=bool(np.all(np.isnan(angles(zero)[0])))
    for name,ok in tests.items():print('PASS' if ok else 'FAIL',name)
    assert all(tests.values()),'Synthetic test failure'
    print('SYNTHETIC',len(tests),'/',len(tests),'PASS')

if __name__=='__main__':verify()
