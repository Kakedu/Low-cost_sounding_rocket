"""ERA5 column fits using the JSBSim hydrostatic temperature-offset atmosphere.

Geopotential heights are converted to geometric heights. Local launch hours
retain the original fixed UTC offsets. Persistence includes midnight pairs.
"""
import csv,datetime,json,sys
from pathlib import Path
import numpy as np,netCDF4
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from atmosphere_model import troposphere_density_kg_m3,EARTH_RADIUS_M
HERE=Path(__file__).resolve().parent
G0=9.80665

def geometric_from_geopotential(H):
    return EARTH_RADIUS_M*np.asarray(H)/(EARTH_RADIUS_M-np.asarray(H))

def fit_dT(z,T,p,pad):
    order=np.argsort(z);z,T,p=np.asarray(z)[order],np.asarray(T)[order],np.asarray(p)[order]
    hs=np.arange(pad,pad+3001,100.0)
    rho=np.exp(np.interp(hs,z,np.log(p)))/(287.05287*np.interp(hs,z,T))
    grid=np.arange(-60,80.001,.5)
    err=np.sum((troposphere_density_kg_m3(hs[None,:],grid[:,None])-rho[None,:])**2,axis=1)
    i=int(np.argmin(err))
    if i in (0,len(grid)-1):raise ValueError('Temperature fit reached search boundary')
    refined=np.arange(grid[i]-.5,grid[i]+.5001,.02)
    err=np.sum((troposphere_density_kg_m3(hs[None,:],refined[:,None])-rho[None,:])**2,axis=1)
    return float(refined[int(np.argmin(err))])

def site(fn,lat,lon,pad,utc_off,hours=(10,17)):
    with netCDF4.Dataset(HERE/fn) as d:
        la,lo=d['latitude'][:],d['longitude'][:]
        lon_q=lon%360 if lo.max()>180 else lon
        i=int(np.argmin(abs(la-lat)));j=int(np.argmin(abs(lo-lon_q)))
        times=netCDF4.num2date(d['time'][:],d['time'].units)
        p=np.asarray(d['level'][:])*100.
        zz=np.asarray(d['z'][:,:,i,j])/G0;tt=np.asarray(d['t'][:,:,i,j])
        out=[]
        for k,t in enumerate(times):
            if not hours[0]<=(t.hour+utc_off)%24<=hours[1]:continue
            utc=datetime.datetime(t.year,t.month,t.day,t.hour,t.minute,t.second)
            out.append((utc,fit_dT(geometric_from_geopotential(zz[k]),tt[k],p,pad)))
    return out

def stats(values):
    v=np.asarray(values,dtype=float)
    if len(v)<2:raise ValueError('Insufficient observations')
    return dict(n=len(v),mean=float(v.mean()),sd=float(v.std(ddof=1)),
                rms_about_isa=float(np.sqrt((v*v).mean())),
                p5=float(np.percentile(v,5)),p95=float(np.percentile(v,95)))

def persistence(series,lag_h):
    mapping=dict(series)
    return stats([x-mapping[t-datetime.timedelta(hours=lag_h)] for t,x in series
                  if t-datetime.timedelta(hours=lag_h) in mapping])

def main():
    eu=site('EuroC_pressure_levels_reanalysis_2001-2021.nc',39.388692,-8.287814,150.,1)
    sp=site('SpaceportAmerica_2018_ERA-5.nc',32.939377,-106.911986,1401.,-6)
    sp_all=site('SpaceportAmerica_2018_ERA-5.nc',32.939377,-106.911986,1401.,-6,hours=(0,23))
    res=dict(euroc_october=stats([x for _,x in eu]),spaceport_2018_all=stats([x for _,x in sp]),
             spaceport_june=stats([x for t,x in sp if t.month==6]),
             pooled_both_sites=stats([x for _,x in eu]+[x for _,x in sp]))
    with netCDF4.Dataset(HERE/'EuroC_single_level_reanalysis_2002_2021.nc') as d:
        la,lo=d['latitude'][:],d['longitude'][:]
        query_lon=-8.287814%360 if lo.max()>180 else -8.287814
        i=int(np.argmin(abs(la-39.388692)));j=int(np.argmin(abs(lo-query_lon)))
        times=netCDF4.num2date(d['time'][:],d['time'].units);t2=np.asarray(d['t2m'][:,i,j])
        mapping={datetime.datetime(t.year,t.month,t.day,t.hour):float(t2[k]) for k,t in enumerate(times)}
    H=150.*EARTH_RADIUS_M/(EARTH_RADIUS_M+150.)
    res['euroc_pad_estimate_error']=stats([mapping[t]-(288.15-.0065*H)-dT for t,dT in eu if t in mapping])
    res['euroc_3h_persistence_error']=persistence(eu,3)
    res['spaceport_6h_persistence_error']=persistence(sp_all,6)
    (HERE/'atmos_stats.json').write_text(json.dumps(res,indent=2),encoding='utf8')
    with (HERE/'atmos_fitted_samples.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['site','utc','delta_T_K','selection'])
        for name,data,selection in [('euroc',eu,'local10_17_fixedUTCplus1'),('spaceport',sp_all,'allhours_fixedUTCminus6')]:
            writer.writerows((name,t.isoformat(),x,selection) for t,x in data)
    (HERE/'atmos_method.json').write_text(json.dumps(dict(
        atmosphere='JSBSim fixed sea-level pressure and temperature-offset hydrostatic troposphere',
        altitude='geometric; ERA5 geopotential converted with radius 6356766m',
        fit='unweighted density least squares from pad to pad+3000m at 100m spacing',
        temperature_search_K=[-60,80],quantile='numpy linear percentile',
        persistence='exact timedifference including midnight',
        site_selection='original fixed local launch windows; site-year subsets are not a universal prior'),indent=2),encoding='utf8')
    for key,value in res.items():print(key,{a:round(b,3) if isinstance(b,float) else b for a,b in value.items()},flush=True)

if __name__=='__main__':main()
