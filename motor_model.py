"""Time-linear RASP curves converted to convergent JSBSim solid-motor tables.

Isp scaling keeps propellant mass and thrust amplitude fixed and stretches
burn duration. start_time is physical elapsed time, including that stretch.
The tiny endpoint regularization tends to zero with the sampling interval.
"""
import bisect
import math
from pathlib import Path

LBF = 4.4482216152605
LB = 0.45359237
# A linear mass ramp from r*F1 to F1 has duration proportional to log(1/r)/(1-r).
# Match the time-linear endpoint ramp duration: log(1/r)=2*(1-r).
# Nontrivial root in (0,1), derived mathematically rather than fitted to apogee.
END_FRACTION = 0.20318786997997995


def read_rasp(path):
    lines = [line.split() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith(';')]
    if len(lines) < 3 or len(lines[0]) < 7:
        raise ValueError(f'Invalid RASP curve: {path}')
    header = lines[0]
    points = [(float(t), float(f)) for t, f in lines[1:]]
    if points[0][0] > 0:
        points.insert(0, (0.0, 0.0))
    validate_points(points)
    prop_kg, total_kg = float(header[4]), float(header[5])
    if not 0 < prop_kg < total_kg:
        raise ValueError('Invalid motor masses')
    return header, points, prop_kg, total_kg


def validate_points(points):
    if points[0] != (0.0, 0.0) or points[-1][1] != 0.0:
        raise ValueError('Curve must start at (0, 0) and end with zero thrust')
    if any(t1 <= t0 or min(f0, f1) < 0 or f0 + f1 == 0
           for (t0, f0), (t1, f1) in zip(points, points[1:])):
        raise ValueError('Times must increase; negative thrust and zero-flow plateaus are unsupported')


def thrust_at(points, t, isp_scale=1.0):
    t = t / isp_scale
    if t <= 0 or t >= points[-1][0]:
        return 0.0
    i = bisect.bisect_right([x[0] for x in points], t) - 1
    t0, f0 = points[i]
    t1, f1 = points[i + 1]
    return f0 + (f1 - f0) * (t - t0) / (t1 - t0)


def impulse_to(points, t):
    """Exact integral of the piecewise time-linear input, through base time t."""
    result = 0.0
    for (t0, f0), (t1, f1) in zip(points, points[1:]):
        if t <= t0:
            break
        end = min(t, t1)
        fend = f0 + (f1 - f0) * (end - t0) / (t1 - t0)
        result += (end - t0) * (f0 + fend) / 2
    return result


def burned_mass(points, prop_kg, t, isp_scale=1.0):
    return prop_kg * impulse_to(points, t / isp_scale) / impulse_to(points, points[-1][0])


def write_motor(points, prop_kg, output, sample_dt=0.002, name='motor',
                start_time=0.0, isp_scale=1.0):
    """Write an engine; tank contents must equal returned remaining_prop_lb.

    points are (seconds, newtons) from ignition. sample_dt is the physical
    resampling interval. A rail-exit model uses start_time=t_exit and starts
    this engine with zero cumulative expended mass, without a buildup delay.
    """
    points = [(float(t), float(f)) for t, f in points]
    validate_points(points)
    if not 0 < sample_dt <= 0.05 or prop_kg <= 0 or isp_scale <= 0:
        raise ValueError('Invalid sampling interval, propellant mass or Isp scale')
    burn_end = points[-1][0] * isp_scale
    if not 0 <= start_time < burn_end:
        raise ValueError('start_time must precede burnout')
    impulse = impulse_to(points, points[-1][0])
    isp = impulse * isp_scale / (LBF * (prop_kg / LB))
    burned_kg = burned_mass(points, prop_kg, start_time, isp_scale)
    remaining_lb = (prop_kg - burned_kg) / LB
    samples = sorted(set(round(t, 10) for t in
        ([t * isp_scale for t, _ in points if t * isp_scale > start_time]
         + [start_time + i * sample_dt for i in range(math.ceil((burn_end-start_time) / sample_dt))]
         + [start_time, burn_end])))
    raw = [(t, thrust_at(points, t, isp_scale)) for t in samples]
    rows = []
    for i, (t, f) in enumerate(raw):
        p = (burned_mass(points, prop_kg, t, isp_scale) - burned_kg) / LB
        effective_f = f
        if i == 0 and start_time == 0:
            effective_f = END_FRACTION * raw[1][1]
        elif i == len(raw)-1:
            effective_f = END_FRACTION * raw[-2][1]
        rows.append((p, effective_f / LBF, t))
    rows[-1] = (remaining_lb, rows[-1][1], rows[-1][2])
    if any(p1 <= p0 for (p0, _, _), (p1, _, _) in zip(rows, rows[1:])):
        raise ValueError('Mass table failed monotonicity')
    body = '\n'.join(f'      {p:.15g} {f:.15g} <!-- source t={t:.9g} s -->' for p, f, t in rows)
    xml = f'''<?xml version="1.0"?>
<!-- {name}: source impulse {impulse:.9g} N.s; sampling {sample_dt:g} s.
     Fixed propellant mass; Isp scale {isp_scale:g}; source start {start_time:g} s.
     Resolution-dependent end regularization; see motor verification. -->
<rocket_engine name="{name}">
  <isp> {isp:.15g} </isp>
  <builduptime> 0 </builduptime>
  <thrust_table name="propulsion/thrust_prop_remain" type="internal">
    <tableData>
{body}
    </tableData>
  </thrust_table>
</rocket_engine>
'''
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(xml, encoding='utf-8')
    return dict(output=str(output.resolve()), impulse_Ns=impulse * isp_scale,
                source_impulse_Ns=impulse, prop_kg=prop_kg, prop_lb=prop_kg/LB,
                burned_kg=burned_kg, remaining_prop_lb=remaining_lb, isp_s=isp,
                isp_scale=isp_scale, start_time_s=start_time, sample_dt_s=sample_dt,
                rows=len(rows), burn_end_s=burn_end,
                initial_regularized_thrust_N=rows[0][1]*LBF,
                terminal_regularized_thrust_N=rows[-1][1]*LBF)


def convert_rasp(path, output, sample_dt=0.002, name='M1670', start_time=0.0, isp_scale=1.0):
    header, points, prop_kg, total_kg = read_rasp(path)
    info = write_motor(points, prop_kg, output, sample_dt, name, start_time, isp_scale)
    info.update(input=str(Path(path).resolve()), total_kg=total_kg,
                hardware_kg=total_kg-prop_kg, rasp_name=header[0])
    return info
