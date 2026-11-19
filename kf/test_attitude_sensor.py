"""Physical sign/projection checks for the measured-attitude SIL model."""
import math
from attitude_sensor import AttitudeSensor, quaternion_from_euler, up_projection

def close(a, b):
    assert abs(a-b) < 1e-9, (a,b)

def main():
    g=9.80665
    # Resting flat: upward accelerometer specific force cancels gravity.
    close(up_projection(quaternion_from_euler(0,0,0),[0,0,-g])-g,0)
    for angle in [70,80,89,90]:
        theta=math.radians(angle); thrust_accel=25.0
        close(up_projection(quaternion_from_euler(0,theta,0),[thrust_accel,0,0])-g,
              thrust_accel*math.sin(theta)-g)
    # Body pitch rate integrates to a 30 deg upward pitch in one second.
    attitude=AttitudeSensor(0,0,0,.02,seed=7,enabled=False)
    for _ in range(50): attitude.step([0,math.radians(30),0])
    close(attitude.specific_force_up([2,0,0]),1.0)
    first=AttitudeSensor(0,0,0,.02,seed=7)
    second=AttitudeSensor(0,0,0,.02,seed=7)
    for _ in range(50):
        first.step([.01,.02,.03]); second.step([.01,.02,.03])
    assert first.q==second.q and first.bias==second.bias
    print("PASS: rest, vertical/inclined force, gyro integration and paired noise")

if __name__=="__main__": main()
