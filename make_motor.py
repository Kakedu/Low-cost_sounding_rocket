"""Convert time-linear RASP thrust to a convergent mass-indexed JSBSim table.

The curve is resampled because a time-linear segment is nonlinear in burned
mass. End values use the resolution-dependent fraction r satisfying
log(1/r)=2(1-r), preserving the small end-interval durations without a
zero-flow equilibrium. The end regularization vanishes as sampling tends
to zero. Tanks must use the returned prop_lb without three-decimal rounding.
"""

import argparse
import json
from pathlib import Path
from motor_model import convert_rasp, read_rasp, thrust_at, LBF, LB


def main():
    ap = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    ap.add_argument('--input', default=str(here / 'M1670.eng'))
    ap.add_argument('--output', default=str(here / 'engine/M1670.xml'))
    ap.add_argument('--sample-dt', type=float, default=0.002)
    ap.add_argument('--name', default='M1670')
    ap.add_argument('--json')
    a = ap.parse_args()
    info = convert_rasp(a.input, a.output, a.sample_dt, a.name)
    print(json.dumps(info, indent=2))
    if a.json:
        Path(a.json).write_text(json.dumps(info, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
