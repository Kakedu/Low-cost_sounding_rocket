"""Rebuild 2022 regional weather from the archived NOAA IGRA excerpt.

Thermodynamic and wind levels have separate validity masks: a missing wind
measurement must never discard a valid pressure/temperature observation.
The whole profiles use their nominal observation times for interpolation;
release times and the regional-proxy limitations remain in the metadata.
"""
from pathlib import Path
import argparse,csv,io,zipfile,json,math,sys,hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from atmosphere_model import EARTH_RADIUS_M,troposphere_density_kg_m3

def records_from_lines(lines):
    wanted={'2022062412','2022062500'}
    records={};raw={};key=None
    for line in lines:
        if line.startswith('#'):
            key=f'{line[13:17]}{line[18:20]}{line[21:23]}{line[24:26]}'
            if key in wanted:
                records[key]={'thermodynamic':[],'wind':[],
                              'release_time_HHMM_UTC':line[27:31]}
                raw[key]=[line]
        elif key in wanted:
            raw[key].append(line)
            p=int(line[9:15]);H=int(line[16:21]);T=int(line[22:27])
            wd=int(line[40:45]);ws=int(line[46:51])
            if H in (-9999,-8888):continue
            h=EARTH_RADIUS_M*H/(EARTH_RADIUS_M-H)
            if p>0 and T not in (-9999,-8888):
                records[key]['thermodynamic'].append((h,T*.1+273.15,p))
            if 0<=wd<=360 and ws>=0:
                east=-ws*.1*math.sin(math.radians(wd))
                north=-ws*.1*math.cos(math.radians(wd))
                records[key]['wind'].append((h,east,north))
    if set(records)!=wanted:raise ValueError('Missing bracketing sounding')
    return records,raw

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--archive',type=Path)
    args=parser.parse_args();data=HERE/'data'
    excerpt=data/'prometheus_2022_igra_excerpt.txt'
    if args.archive:
        with zipfile.ZipFile(args.archive) as z,io.TextIOWrapper(z.open(z.namelist()[0]),encoding='ascii') as f:
            records,raw=records_from_lines(f)
        excerpt.write_text(''.join(''.join(raw[k]) for k in sorted(raw)),encoding='ascii')
    else:
        with excerpt.open(encoding='ascii') as f:records,raw=records_from_lines(f)
    hs=np.arange(1401.,6201.,100.);profiles=[]
    for key in sorted(records):
        therm=np.array(sorted(records[key]['thermodynamic']))
        wind=np.array(sorted(records[key]['wind']))
        for rows in (therm,wind):
            if len(rows)<2 or rows[0,0]>=hs[0] or rows[-1,0]<=hs[-1]:
                raise ValueError('Sounding does not bracket fitting column')
        profile=np.column_stack((hs,np.interp(hs,wind[:,0],wind[:,1]),
                                np.interp(hs,wind[:,0],wind[:,2]),
                                np.interp(hs,therm[:,0],therm[:,1]),
                                np.exp(np.interp(hs,therm[:,0],np.log(therm[:,2])))))
        profiles.append(profile)
        with (data/f'prometheus_sounding_{key}.csv').open('w',newline='') as f:
            writer=csv.writer(f)
            writer.writerow(['geometric_h_m','wind_east_ms','wind_north_ms','temperature_K','pressure_Pa'])
            writer.writerows(profile)
    weight=(3+17/60)/12
    profile=(1-weight)*profiles[0]+weight*profiles[1]
    profile=np.column_stack((profile,profile[:,4]/(287.05287*profile[:,3])))
    grid=np.arange(-30.,60.001,.01)
    err=np.sum((troposphere_density_kg_m3(hs[:,None],grid[None,:])-profile[:,5,None])**2,axis=0)
    dT=float(grid[int(np.argmin(err))])
    la1,la2=map(math.radians,(31.8728,32.939377));dlo=math.radians(-106.911986+106.6980)
    distance=6371000*2*math.asin(math.sqrt(math.sin((la2-la1)/2)**2+math.cos(la1)*math.cos(la2)*math.sin(dlo/2)**2))
    base=json.loads((HERE/'prometheus_inputs.json').read_text())
    base.update(source='NOAA IGRA2.2 USM00072364 independent 2022 regional profiles; published geometry retained',
        weather_time_utc='2022-06-24 15:17:00',profile=profile.tolist(),
        weather_metadata=dict(station='USM00072364',station_latitude=31.8728,station_longitude=-106.6980,
            station_pad_distance_m=distance,bracketing_times_utc=['2022-06-24 12:00','2022-06-25 00:00'],
            nominal_hour_weight=weight,fit_dT_K=dT,
            source_url='https://www.ncei.noaa.gov/pub/data/igra/data/data-por/USM00072364-data.txt.zip',
            excerpt_sha256=hashlib.sha256(excerpt.read_bytes()).hexdigest(),
            profiles={key:dict(thermodynamic_levels=len(r['thermodynamic']),wind_levels=len(r['wind']),
                        release_time_HHMM_UTC=r['release_time_HHMM_UTC']) for key,r in records.items()},
            processing='Separate valid thermodynamic/wind masks; logarithmic vertical pressure interpolation; meteorological FROM-wind converted to east/north; nominal-hour interpolation',
            limitation='Regional proxy about 120km away, not launch-site weather; whole-profile nominal-time interpolation ignores balloon ascent time; spatial and profile-shape representativeness unquantified'))
    (HERE/'prometheus_inputs_2022.json').write_text(json.dumps(base,indent=2),encoding='utf-8')
    print(json.dumps(base['weather_metadata'],indent=2))

if __name__=='__main__':main()
