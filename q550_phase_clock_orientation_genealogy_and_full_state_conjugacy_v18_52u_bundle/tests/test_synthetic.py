from pathlib import Path
import sys,numpy as np,math
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from energy_math import Params,sigmoid
from clock_math import field_plus,field_external_minus,relabelled_field,FLIP_THETA,FLIP_ALL_PHASES,general_conjugacy_residual,normres
p=Params()
rng=np.random.default_rng(18525220)
match=[];mismatch=[];clock=[];reflect=[]
for i in range(200):
 z=np.zeros(11)
 z[[0,4]]=rng.uniform(-1,1,2);z[[1,3,5,7,9]]=rng.uniform(-.5,.5,5)
 z[[2,6,10]]=rng.uniform(-3,3,3);z[8]=rng.uniform(-2,2)
 k=[1./6,.1665996104][i%2]
 plus=field_plus(z,k,p);minus=field_external_minus(z,k,p)
 rec=relabelled_field(z,k,p)
 exact=FLIP_THETA*field_plus(FLIP_THETA*z,k,p)
 match.append(normres(rec,exact));clock.append(abs(plus[10]+minus[10]))
 mismatch.append(normres(rec,minus))
 reflect.append(normres(general_conjugacy_residual(z,k,p,FLIP_ALL_PHASES)[1],minus))
 assert abs(plus[10]-float(sigmoid(z[8])))<1e-14
 assert abs(minus[10]+float(sigmoid(z[8])))<1e-14
assert max(match)<1e-15 and max(clock)<1e-15
assert max(mismatch)>1e-4 and max(reflect)>1e-4
print('SYNTHETIC_PASS',{'coordinate_relabel_max_relative':max(match),'naive_external_minus_mismatch_max':max(mismatch),'all_angle_flip_mismatch_max':max(reflect)})
