"""Shared fail-fast checks; drivers must also prove their requested endpoint."""
import math

def require(ok, action):
    if not ok:
        raise RuntimeError(f'JSBSim {action} failed')

def run_checked(fdm):
    require(fdm.run(), 'integration')
    if not math.isfinite(fdm['position/h-agl-ft']):
        raise RuntimeError('Non-finite altitude')
