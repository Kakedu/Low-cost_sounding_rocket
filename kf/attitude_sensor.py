"""Attitude estimate from a launch inclinometer and three noisy body gyros.

The numerical spreads are engineering assumptions, not flight measurements:
initial Euler angles 0.1 deg SD each, fixed gyro residual bias 0.02 deg/s SD,
and each 50 Hz gyro sample 0.1 deg/s SD. Rates are body rates in rad/s.
"""
import hashlib
import math
import random

INITIAL_SD_DEG = 0.1
GYRO_BIAS_SD_DEG_S = 0.02
GYRO_SAMPLE_SD_DEG_S = 0.1

def named_rng(seed, name):
    digest = hashlib.sha256(f"{seed}:{name}".encode()).digest()
    return random.Random(int.from_bytes(digest[:8], "little"))

def quaternion_from_euler(phi, theta, psi):
    cp, sp = math.cos(phi / 2), math.sin(phi / 2)
    ct, st = math.cos(theta / 2), math.sin(theta / 2)
    cy, sy = math.cos(psi / 2), math.sin(psi / 2)
    return [cp*ct*cy + sp*st*sy, sp*ct*cy - cp*st*sy,
            cp*st*cy + sp*ct*sy, cp*ct*sy - sp*st*cy]

def up_projection(q, vector):
    w, x, y, z = q
    # Negative row 3 of body-to-NED rotation matrix.
    up = [2*(w*y-x*z), -2*(y*z+w*x), -(1-2*(x*x+y*y))]
    return sum(a*b for a,b in zip(up, vector))

class AttitudeSensor:
    def __init__(self, phi, theta, psi, dt, seed, enabled=True):
        self.dt = dt
        self.rng = named_rng(seed, "gyro_samples")
        initial = named_rng(seed, "initial_attitude")
        bias = named_rng(seed, "gyro_bias")
        angles = [a + math.radians(initial.gauss(0, INITIAL_SD_DEG)) * enabled
                  for a in [phi, theta, psi]]
        self.bias = [math.radians(bias.gauss(0, GYRO_BIAS_SD_DEG_S))*enabled for _ in range(3)]
        self.q = quaternion_from_euler(*angles)
        self.enabled = enabled

    def step(self, body_rates):
        rates = [r + b + math.radians(self.rng.gauss(0, GYRO_SAMPLE_SD_DEG_S))*self.enabled
                 for r,b in zip(body_rates,self.bias)]
        norm = math.sqrt(sum(r*r for r in rates))
        half = norm * self.dt / 2
        scale = math.sin(half)/norm if norm > 1e-12 else self.dt/2
        a,b,c,d = [math.cos(half), *[scale*r for r in rates]]
        w,x,y,z = self.q
        q = [w*a-x*b-y*c-z*d, w*b+x*a+y*d-z*c,
             w*c-x*d+y*a+z*b, w*d+x*c-y*b+z*a]
        mag = math.sqrt(sum(v*v for v in q)); self.q = [v/mag for v in q]

    def specific_force_up(self, body_specific_force):
        return up_projection(self.q, body_specific_force)
