"""Derive the small runtime curve data from the supplied S1, without rounding."""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1]/'literature/manual_inbox/contact_real_data/static/Data Sheet 1.CSV'


def derive():
    data = {}
    with SOURCE.open(encoding='utf-8-sig', newline='') as stream:
        for row in csv.DictReader(stream):
            if row['id'].startswith('P') or row['id'] == 'rangeMean':
                degree = int(row['fitModel'][4:])
                data.setdefault(row['position'], {})[row['id']] = dict(
                    coeff=[float(row[f'p{i}']) for i in range(1, degree+2)],
                    fitMaxDispl=float(row['fitMaxDispl']))
    return data


if __name__ == '__main__':
    (HERE/'remus_s1.json').write_text(json.dumps(derive(), indent=2, allow_nan=False)+'\n', encoding='utf-8')
