"""Pressure-altitude observation model for the JSBSim validation atmosphere.

The log's proprietary conversion is not documented. A standard dry ICAO
tropospheric pressure-altitude conversion, zeroed at the pad, is assumed.
Both pressures come from the same geometric-altitude JSBSim hydrostatic
temperature-offset model used by the flight simulation.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from atmosphere_model import troposphere_pressure_pa,T0_K,P0_PA,LAPSE_K_PER_M,R_AIR,G0

def standard_pressure_altitude_m(pressure_pa):
    return T0_K/LAPSE_K_PER_M*(1-(pressure_pa/P0_PA)**(R_AIR*LAPSE_K_PER_M/G0))

def baro_reading(h_true,dT,pad=150.0):
    p0=troposphere_pressure_pa(pad,dT)
    p1=troposphere_pressure_pa(pad+h_true,dT)
    return float(standard_pressure_altitude_m(p1)-standard_pressure_altitude_m(p0))
