"""Dry troposphere consistent with JSBSim delta-T at fixed sea-level pressure.

Altitude is geometric metres ASL; delta-T is kelvin. This model covers the
troposphere only, as used by the present rocket studies. Scalar or NumPy inputs
are supported. The scientific 101325 Pa sea-level pressure differs from
JSBSim's 2116.228 psf by about 5 ppm.
"""
G0 = 9.80665
R_AIR = 287.05287
EARTH_RADIUS_M = 6356766.0
LAPSE_K_PER_M = 0.0065
P0_PA = 101325.0
T0_K = 288.15

def troposphere_pressure_pa(h_m, dT_K=0.0):
    H = h_m * EARTH_RADIUS_M / (EARTH_RADIUS_M + h_m)
    T0 = T0_K + dT_K
    T = T0 - LAPSE_K_PER_M * H
    return P0_PA * (T / T0) ** (G0 / (R_AIR * LAPSE_K_PER_M))

def troposphere_density_kg_m3(h_m, dT_K=0.0):
    H = h_m * EARTH_RADIUS_M / (EARTH_RADIUS_M + h_m)
    T = T0_K + dT_K - LAPSE_K_PER_M * H
    return troposphere_pressure_pa(h_m, dT_K) / (R_AIR * T)
