"""Flight-informed pressure/GPS column calibration, not independent weather.

GPS heights are geometric above the first fix. The assumed pad geometric ASL
altitude is 1401m. Pressure is averaged over the first five records, while
apex pressure is the minimum in +/-1s of the fitted GPS apex. Timing/port
offsets are not fitted. The diagnostic cannot identify a unique temperature
profile or independently validate the apogee used to derive it.
"""
from pathlib import Path
import csv,hashlib
import numpy as np
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from atmosphere_model import EARTH_RADIUS_M,R_AIR,G0,LAPSE_K_PER_M,T0_K

LOGS={'TeleMetrum':'2022-06-24-serial-5115-flight-0001-TeleMetrum.csv',
      'TeleMega':'2022-06-24-serial-6583-flight-0003-TeleMega.csv'}
PAD_M=1401.

def column_temperature_K(p_pad,p_apex,geometric_dz_m,pad_geometric_ASL_m=PAD_M):
    H=lambda h:EARTH_RADIUS_M*h/(EARTH_RADIUS_M+h)
    hp=H(pad_geometric_ASL_m);ha=H(pad_geometric_ASL_m+geometric_dz_m)
    x=(p_apex/p_pad)**(R_AIR*LAPSE_K_PER_M/G0)
    if not 0<x<1:raise ValueError('Apex pressure must be below positive pad pressure')
    return LAPSE_K_PER_M*(ha-hp)/(1-x)-(T0_K-LAPSE_K_PER_M*hp)

def legacy_column_temperature_K(p_pad,p_apex,geometric_dz_m,pad_geometric_ASL_m=PAD_M):
    """Retain the superseded flat-height diagnostic for the correction log."""
    x=(p_apex/p_pad)**(287.053*.0065/9.80665)
    return .0065*geometric_dz_m/(1-x)-(288.15-.0065*pad_geometric_ASL_m)

def infer_from_logs(data_dir=None):
    data=Path(data_dir) if data_dir is not None else Path(__file__).resolve().parent/'data'
    records=[]
    for name,fn in LOGS.items():
        path=data/fn;rows=list(csv.reader(path.open()));head=[x.strip() for x in rows[0]]
        it=head.index('time');ip=head.index('pressure');ig=max(i for i,x in enumerate(head) if x=='altitude')
        raw=[(float(r[it]),float(r[ip]),float(r[ig])) for r in rows[1:] if len(r)==len(head)]
        ppad=sum(r[1] for r in raw[:5])/5;padgps=raw[0][2];fixes=[];previous=None
        for t,p,g in raw:
            if g!=previous:fixes.append((t,g-padgps));previous=g
        tm=max((r for r in fixes if r[0]>10),key=lambda r:r[1])[0]
        pts=[r for r in fixes if abs(r[0]-tm)<=4];co=np.polyfit([t-tm for t,g in pts],[g for t,g in pts],2)
        if co[0]>=0:raise ValueError('GPS apex fit does not open downwards')
        tv=-co[1]/(2*co[0]);t_apex=tm+tv;g_apex=float(np.polyval(co,tv))
        papex=min(p for t,p,g in raw if abs(t-t_apex)<1)
        records.append(dict(name=name,source_file=fn,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),gps_apex_geometric_m=g_apex,gps_apex_time_s=float(t_apex),gps_fit_points=len(pts),pad_pressure_Pa=ppad,apex_pressure_Pa=papex,pad_gps_first_fix_m=padgps,assumed_pad_geometric_ASL_m=PAD_M,column_dT_K=float(column_temperature_K(ppad,papex,g_apex)),legacy_flat_column_dT_K=float(legacy_column_temperature_K(ppad,papex,g_apex))))
    return dict(method='Geometric-to-geopotential pressure-ratio inversion; mean of the two GPS-apex diagnostics',role='flight-informed calibration',records=records,column_dT_K=sum(r['column_dT_K'] for r in records)/len(records),legacy_flat_column_dT_K=sum(r['legacy_flat_column_dT_K'] for r in records)/len(records))
