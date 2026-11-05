"""Reads the four Neumann et al. 2024 spinal cord indentation archives (Zenodo 14536174) into arrays, nothing else.

Chiaro nanoindenter, sphere about 100 um in radius, 1 kHz. Columns: time s, load uN, indentation nm, cantilever nm, piezo nm.
Files named NOK by the instrument failed its automatic Hertz fit, not necessarily the measurement, so each file is judged from its own curve.
"""
import re
import zipfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RAW = HERE.parents[1]/'literature'/'manual_inbox'/'contact_real_data'
ARCHIVES = dict(rate='Figure 7 & 8 - Loading rate sensitivity test.zip', relax='Figure 12 - Stress relaxation test.zip',
                cyclic='Figure 13 - Cyclic loading test.zip', dma='Figure 14 - Dynamic mechanical analysis.zip')


def _floats(text):
    return [float(x) for x in re.findall(r'-?\d+\.?\d*', text)]


def read(archive, member):
    with zipfile.ZipFile(RAW/archive) as z, z.open(member) as f:
        raw = f.read().decode('latin1')
    head, body = raw.split('Time (s)', 1)
    body = body.split('\n', 1)[1]
    data = np.fromstring(body.replace('\r', ' ').replace('\n', ' ').replace('\t', ' '), sep=' ')
    data = data.reshape(-1, 6)
    meta = dict(k=float(re.search(r'k \(N/m\)\s+([\d.]+)', head)[1]), R=float(re.search(r'Tip radius \(um\)\s+([\d.]+)', head)[1]),
                status=re.search(r'Status\s+(\w+)', head)[1], specimen=re.search(r'\n(\S+)\s*\r?\n', head)[1], name=member.split('/')[-1][:-4],
                tissue='grey' if 'grey' in member else 'white',
                start=_floats(re.search(r'Step absolute start times \(s\)(.*)', head)[1]), end=_floats(re.search(r'Step absolute end times \(s\)(.*)', head)[1]))
    m = re.search(r'Depth \(nm\) ([\d.]+)', head)
    if m:
        meta['depth'] = float(m[1])*1e-3
        meta['freq'] = [float(x) for x in re.findall(r'Freq ([\d.]+)', head)]
        meta['amp'] = [float(x)*1e-3 for x in re.findall(r'Amp ([\d.]+)', head)]
        meta['dma_start'] = _floats(re.search(r'DMA absolute start times \(s\)(.*)', head)[1])
        meta['dma_end'] = _floats(re.search(r'DMA absolute end times \(s\)(.*)', head)[1])
    else:
        meta['profile'] = [(float(a)*1e-3, float(b)) for a, b in re.findall(r'D\[Z\d+\] \(nm\)\s+([\d.]+)\s+t\[\d+\] \(s\)\s+([\d.]+)', head)]
    return meta, dict(t=data[:, 0], load=data[:, 1], ind=data[:, 2]*1e-3, cant=data[:, 3]*1e-3, piezo=data[:, 4]*1e-3)


def files(kind):
    with zipfile.ZipFile(RAW/ARCHIVES[kind]) as z:
        return ARCHIVES[kind], sorted(n for n in z.namelist() if n.endswith('.txt'))
