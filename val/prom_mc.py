"""Run the corrected prometheus uncertainty study; full suite: validation_mc.py."""
import argparse
from validation_mc import monte_carlo

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('n', nargs='?', type=int, default=200)
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    if args.n < 2 or args.workers < 1:
        parser.error('n must be at least 2 and workers must be positive')
    result = monte_carlo('prometheus', args.n, args.workers)
    print({k: v for k, v in result.items() if k not in ('samples', 'assumptions')})
