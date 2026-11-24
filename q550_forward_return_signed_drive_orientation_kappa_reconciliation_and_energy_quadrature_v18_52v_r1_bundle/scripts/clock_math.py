"""Physical-time Q550 clock maps. No ODE integration or fitted sign."""
import numpy as np
from phase_energy import physical_rhs_with_signed_external_phase

FLIP_THETA=np.array([1,1,1,1,1,1,1,1,1,1,-1],float)
FLIP_ALL_PHASES=np.array([1,1,-1,-1,1,1,-1,-1,1,1,-1],float)

def field_plus(y,k,p):
    return physical_rhs_with_signed_external_phase(y,k,p,float(y[10]),sigma=+1)

def field_external_minus(y,k,p):
    return physical_rhs_with_signed_external_phase(y,k,p,float(y[10]),sigma=-1)

def relabelled_field(y_ext,k,p):
    """Coordinate relabel Theta_ext=-Theta_engine; delta becomes sum+Theta_ext."""
    eng=FLIP_THETA*np.asarray(y_ext,float)
    return FLIP_THETA*field_plus(eng,k,p)

def general_conjugacy_residual(y,k,p,M):
    """Compare the archived minus-clock model with pushforward of original plus-clock field."""
    x=M*np.asarray(y,float)
    push=M*field_plus(x,k,p)
    observed=field_external_minus(y,k,p)
    return push-observed,push,observed

def angular_accelerations(Y,k,p,phase,theta_sign=-1):
    """Predicted omega derivatives using either sign with recorded external phase."""
    A=np.asarray(Y,float)
    r1=np.exp(A[:,0]);r2=np.exp(A[:,4]);x=1./(1.+np.exp(-A[:,8]))
    D=A[:,2]+A[:,6]-p.phi0+theta_sign*np.asarray(phase,float)
    o1=k*x*(r2/r1)*np.sin(D)-2*A[:,1]*A[:,3]
    o2=k*x*(r1/r2)*np.sin(D)-2*A[:,5]*A[:,7]
    return np.stack([o1,o2],axis=1)

def normres(a,b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/max(1.,np.linalg.norm(np.asarray(a)),np.linalg.norm(np.asarray(b))))
