"""Regenerate the current Prometheus geometric/barometric comparison plot."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).with_name('gps_check.py')), run_name='__main__')
