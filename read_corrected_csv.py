"""Read legacy CSVs with corrected pitch-rate field names; values stay rad/s."""
import csv
def read_corrected_csv(path):
    with open(path, newline='', encoding='utf-8') as stream:
        reader=csv.reader(stream)
        names=[x.strip().replace('Pitch rate (deg/s)', 'Pitch rate (rad/s)') for x in next(reader)]
        yield from (dict(zip(names,row)) for row in reader)
