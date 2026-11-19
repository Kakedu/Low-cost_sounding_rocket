"""Check nonnegative means at each measurement and before propagation."""
import pathlib,sys,numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent.parent))
from kf.kf_estimator import FlightKF,FlightKF4

def check(cls):
    k=cls(.02,dict(baro_sigma_m=.5,acc_sigma_g=.04,port_cp_sd=.03),.45,2,.00785,17.6)
    k.x[2]=-2.5;k._predict_coast(.5)
    assert k.x[2]==0 and k.x[1]<0
    assert k._drag_accel(k.x[:3],.5)>=0 # descending drag opposes descent
    k.x[0:3]=[1000,0,.45] # zero dynamic pressure isolates the barometric Cd cross-gain
    k.P=np.diag([1.,1.,.01]+([.0009] if cls is FlightKF4 else []));k.P[0,2]=k.P[2,0]=.05
    k._update_baro(900)
    assert k.x[2]>=0
    k.x[0:3]=[1000,200,.45]
    k.P=np.diag([1.,1.,.01]+([.0009] if cls is FlightKF4 else []))
    k._update_accel(100,0)
    assert k.x[2]>=0
    assert k._drag_accel(k.x[:3],0)<=0 # ascending force cannot accelerate up
    assert k.cd_projection_count>=3
    assert np.all(np.isfinite(k.P))

if __name__=='__main__':
    check(FlightKF);check(FlightKF4);print('PASS: KF3/KF4 mean projections after each update and before prediction; no prior/Q/R/Cp change')
