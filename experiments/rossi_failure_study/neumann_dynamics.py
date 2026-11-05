"""Identifies and cross-validates dynamic contact mechanisms on the four Neumann et al. 2024 protocols: loading rate, relaxation, cyclic, DMA.

Forward model, per file, with the measured displacement as input:
    F = S * ( w_e u + sum w_i q_i + w_v n d^(n-1) d' + w_hc u d' ),  u = d^n,  q_i' = u' - q_i/tau_i

The shape w is shared within a tissue; only S is free per file, because each file is another site or specimen. Contact is judged on the measured force.
Files that the instrument marked NOK failed its own Hertz fit, so the contact point is refitted here for every ramp from F^(2/3) against displacement.
"""
import json
import math
import os
import pickle
import tempfile
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, nnls

import neumann_data as N

HERE = Path(__file__).resolve().parent
TAU = np.geomspace(1.e-3, 1.e3, 19)
EXP = 1.6        # elastic exponent: 1.5 apparent in the ramps plus the slope of memory, and .55-.65+1 from stiffness against depth in the DMA
APPARENT = 1.5
CACHE = Path(os.environ.get('NEUMANN_CACHE', tempfile.gettempdir()))/'neumann_prepared.pkl'


def contact_point(z, f):
    """Fits f^(1/1.5) linear in z on the first loading, between 20 and 80 percent of the peak; returns the contact position and the fit error."""
    pk = np.argmax(z)
    z1, f1 = z[:pk+1], f[:pk+1]
    fmax = f1.max()
    i_hi = np.argmax(f1 >= .8*fmax)
    i_lo = np.where(f1[:i_hi] <= .2*fmax)[0]
    i_lo = i_lo[-1] if len(i_lo) else 0
    y = np.maximum(f1[i_lo:i_hi+1], 0.)**(1./APPARENT)
    a, b = np.polyfit(z1[i_lo:i_hi+1], y, 1)
    return -b/a, float(np.std(y-(a*z1[i_lo:i_hi+1]+b))/(y.max()+1e-12))


def exponent(z, f, z0):
    """Local exponent of force against depth on the first loading above 30 percent of the peak depth."""
    pk = np.argmax(z); d = z[:pk+1]-z0; ff = f[:pk+1]
    q = (d > .3*d[-1]) & (ff > .05*ff.max())
    return float(np.polyfit(np.log(d[q]), np.log(ff[q]), 1)[0]) if q.sum() > 5 else float('nan')


def resample(t, z, f, events):
    """Averages into bins, fine after each event in the protocol and coarse when nothing happens, so a 600 s hold does not outweigh the ramp."""
    edges = [t[0]]
    ev = np.sort(np.asarray(events)); ev = ev[(ev > t[0]) & (ev < t[-1])]
    while edges[-1] < t[-1]:
        now = edges[-1]
        past = ev[ev <= now]; next_ev = ev[ev > now]
        since = now-past[-1] if len(past) else 1e9
        h = float(np.clip(.03*since, .005, 8.))
        if len(next_ev):
            h = min(h, max(next_ev[0]-now, .005))
        edges.append(now+h)
    idx = np.searchsorted(edges, t, side='right')-1
    n = idx.max()+1
    count = np.bincount(idx, minlength=n).astype(float); keep = count > 0
    mean = lambda x: (np.bincount(idx, x, minlength=n)[keep]/count[keep])
    return mean(t), mean(z), mean(f)


def prepare_ramp(kind, member):
    arc, _ = N.files(kind)
    m, d = N.read(arc, member)
    t, f, z = d['t'], d['load'], d['piezo']-d['cant']
    f = f-np.median(f[:200])       # free liquid at the start of the record
    prof = m['profile']
    pk = int(np.argmax(z[:int(len(z)*.8)]))
    if kind == 'cyclic':       # the first loading alone, so the record up to the first turn
        first = np.argmax(np.diff(z, prepend=z[0])[500:] < -5.e-3)+500
        pk = int(np.argmax(z[:first+1]))
    z0, err = contact_point(z[:pk+1], f[:pk+1])
    tc = t[np.argmax(z >= z0)]
    if kind == 'relax':
        drop = pk+np.argmax(z[pk:] < z[pk:pk+2000].mean()-1.)
    elif kind == 'rate':
        drop = pk+np.argmax(z[pk:] < z[pk:pk+500].mean()-1.)
    else:
        drop = pk
    events = [tc, t[pk], t[drop]]
    if kind == 'cyclic':
        events = [tc]+[t[pk]+i*5.018 for i in range(10)]
    tb, zb, fb = resample(t, z, f, events)
    return dict(kind=kind, name=m['name'], tissue=m['tissue'], specimen=m['specimen'], status=m['status'], k=m['k'], R=m['R'], t=tb, d=zb-z0, f=fb,
                fit_error=err, n_local=exponent(z[:pk+1], f[:pk+1], z0), speed=float(prof[0][0]/prof[0][1]), peak_force=float(f.max()), depth_peak=float(z[pk]-z0), t_peak=float(t[pk]),
                t_drop=float(t[drop]), t_contact=float(tc))


def prepare_dma(member):
    arc, _ = N.files('dma')
    m, d = N.read(arc, member)
    rows = []
    for f, s, e in zip(m['freq'], m['dma_start'], m['dma_end']):
        q = (d['t'] >= s+1./f) & (d['t'] <= e)
        t = d['t'][q]-d['t'][q].mean()
        a = np.c_[np.sin(2*np.pi*f*t), np.cos(2*np.pi*f*t), np.ones_like(t), t, t*t]
        cf = np.linalg.lstsq(a, d['load'][q], rcond=None)[0]; cd = np.linalg.lstsq(a, d['ind'][q], rcond=None)[0]
        ff = complex(cf[0], cf[1]); dd = complex(cd[0], cd[1])
        rows.append(dict(f=f, K=ff/dd, F=abs(ff), D=abs(dd), d0=float(cd[2]), F0=float(cf[2]), snr=abs(ff)/float((d['load'][q]-a@cf).std())))
    return dict(kind='dma', name=m['name'], tissue=m['tissue'], specimen=m['specimen'], status=m['status'], k=m['k'], R=m['R'], target=m['depth'], rows=rows)


def prepare():
    if CACHE.exists():
        return pickle.loads(CACHE.read_bytes())
    jobs = [(k, n) for k in ('rate', 'relax', 'cyclic') for n in N.files(k)[1]]
    with ProcessPoolExecutor(max_workers=12) as ex:
        data = list(ex.map(prepare_ramp, [j[0] for j in jobs], [j[1] for j in jobs], chunksize=1))
        data += list(ex.map(prepare_dma, N.files('dma')[1], chunksize=1))
    CACHE.write_bytes(pickle.dumps(data))
    return data


def basis(r):
    """Columns for one time record: e, the branches on the grid, stiffness-proportional viscous n d^(n-1) d' (weight in s), Hunt-Crossley u d' (weight lambda in s/um)."""
    t, d = r['t'], np.maximum(r['d'], 0.)
    u = d**EXP
    h = np.diff(t, prepend=t[0]); h[0] = h[1]
    du = np.diff(u, prepend=u[0]); dd = np.gradient(r['d'], t)
    cols = [u]
    for tau in TAU:
        q = np.zeros_like(u); e = np.exp(-h/tau); g = tau/h*(1.-e)
        for k in range(1, len(u)):
            q[k] = e[k]*q[k-1]+du[k]*g[k]
        cols.append(q)
    cols.append(EXP*np.maximum(d, 1e-9)**(EXP-1.)*dd)
    cols.append(u*dd)
    return np.array(cols).T


def in_contact(r):
    """Depth positive and force clearly positive: once the surface has come off the tip, or sticks to it, none of these models describe the force."""
    return (r['d'] > 0.) & (r['f'] > .02*r['peak_force'])


def dma_row(r, row):
    d0 = row['d0'] if r['name'] != 'PSC1501 Indentation_032' else r['target']       # this file has a wrong contact point, by its load and its target
    w = 2.*np.pi*row['f']
    grow = EXP*d0**(EXP-1.)
    return np.r_[grow, grow*(1j*w*TAU/(1.+1j*w*TAU)), 1j*w*grow, 1j*w*d0**EXP]


class Data:
    def __init__(self, records):
        self.t = [r for r in records if r['kind'] != 'dma']
        self.m = [r for r in records if r['kind'] == 'dma']
        self.X = [basis(r) for r in self.t]
        self.mask = [in_contact(r) for r in self.t]
        self.D = [np.array([dma_row(r, x) for x in r['rows']]) for r in self.m]
        self.K = [np.array([x['K'] for x in r['rows']]) for r in self.m]


def objective(data, w, protocols):
    """One number for every family: mean over protocols of the squared relative error, each file one vote, each file with its own scale."""
    total = 0.
    for proto in protocols:
        parts = []
        if proto == 'dma':
            for Dm, K in zip(data.D, data.K):
                m = (Dm@w)/np.abs(K); y = K/np.abs(K)
                s = float(np.real(np.vdot(m, y)))/max(float(np.real(np.vdot(m, m))), 1e-30)
                parts.append(float(np.mean(np.abs(s*m-y)**2)))
        else:
            for r, X, c in zip(data.t, data.X, data.mask):
                if r['kind'] != proto or c.sum() < 5:
                    continue
                x = X[c]@w; f = r['f'][c]
                s = float(f@x)/max(float(x@x), 1e-30)
                parts.append(float(np.mean((s*x-f)**2))/r['peak_force']**2)
        total += float(np.mean(parts))
    return total/len(protocols)


E_COL, FIRST, KV, HC = 0, 1, 1+len(TAU), 2+len(TAU)


def alternate(data, protocols, w, active, iters=40):
    """Same objective for a shape free on the active columns: scale per file, then nonnegative least squares, keeping the best seen."""
    active = np.asarray(active)
    best, best_j = w.copy(), objective(data, w, protocols)
    for _ in range(iters):
        A, b = [], []
        for proto in protocols:
            nf = sum(1 for r in data.t+data.m if r['kind'] == proto)
            if proto == 'dma':
                for Dm, K in zip(data.D, data.K):
                    m = (Dm@w)/np.abs(K); y = K/np.abs(K)
                    s = float(np.real(np.vdot(m, y)))/max(float(np.real(np.vdot(m, m))), 1e-30)
                    wt = 1./math.sqrt(nf*len(K))
                    M = s*(Dm.T/np.abs(K)).T*wt
                    A += [M.real, M.imag]; b += [y.real*wt, y.imag*wt]
            else:
                for r, X, c in zip(data.t, data.X, data.mask):
                    if r['kind'] != proto or c.sum() < 5:
                        continue
                    x = X[c]; f = r['f'][c]; xw = x@w
                    s = float(f@xw)/max(float(xw@xw), 1e-30)
                    wt = 1./(r['peak_force']*math.sqrt(nf*len(f)))
                    A.append(s*wt*x); b.append(wt*f)
        A = np.vstack(A); b = np.concatenate(b)
        sol, _ = nnls(A[:, active], b, maxiter=20000)
        w = np.zeros_like(w); w[active] = sol
        w = w/max(w[:KV].sum(), 1e-12)
        j = objective(data, w, protocols)
        if j < best_j-1e-9:
            best, best_j = w.copy(), j
    return best


def score(data, w, protocols=('rate', 'relax', 'cyclic', 'dma')):
    """Prediction errors with only S fitted to each file. Time records: error against the file's peak, the prediction clipped at no force."""
    out = {}
    for proto in protocols:
        if proto == 'dma':
            ph, gain = [], []
            for Dm, K in zip(data.D, data.K):
                m = Dm@w
                s = float(np.exp(np.mean(np.log(np.abs(K))-np.log(np.abs(m)+1e-30))))
                ph.append(np.degrees(np.angle(m))-np.degrees(np.angle(K)))
                gain.append(np.log(np.abs(m)*s)-np.log(np.abs(K)))
            ph, gain = np.array(ph), np.array(gain)
            out['dma_phase_rmse_deg'] = float(np.sqrt(np.mean(ph**2)))
            out['dma_phase_error_20Hz_deg'] = float(np.mean(ph[:, -1]))
            out['dma_phase_error_05Hz_deg'] = float(np.mean(ph[:, 0]))
            out['dma_gain_rmse_pct'] = float(100.*np.sqrt(np.mean(gain**2)))
            continue
        err = []
        for r, X, c in zip(data.t, data.X, data.mask):
            if r['kind'] != proto or c.sum() < 5:
                continue
            x = X[c]@w; f = r['f'][c]
            s = float(f@x)/max(float(x@x), 1e-30)
            err.append(float(np.sqrt(np.mean((np.maximum(s*x, 0.)-f)**2))/r['peak_force']))
        out[proto+'_nrmse_pct'] = float(100.*np.mean(err)) if err else float('nan')
    return out


def zero():
    return np.zeros(HC+1)


def power_shape(beta, fe, viscous=0.):
    w = zero(); a = TAU**(-beta); w[FIRST:KV] = a/a.sum()*(1.-fe); w[E_COL] = fe; w[KV] = viscous
    return w


def search(data, protocols, make, grids):
    """Grid, then simplex on the same objective."""
    best = None
    for point in product(*grids):
        j = objective(data, make(*point), protocols)
        if best is None or j < best[0]:
            best = (j, point)
    f = lambda x: objective(data, make(*[max(v, 0.) for v in x]), protocols)
    res = minimize(f, np.array(best[1])+1e-4, method='Nelder-Mead', options=dict(maxiter=300, xatol=1e-4, fatol=1e-7))
    return make(*[max(v, 0.) for v in res.x])


def nearest(x):
    return int(np.argmin(np.abs(np.log(TAU/x))))


def free_on(data, protocols, active, start=None):
    active = list(active)
    w = power_shape(.15, .3) if start is None else start
    keep = zero(); keep[active] = w[active]
    if not np.any(keep > 0.):
        keep[active] = 1./len(active)
    return alternate(data, protocols, keep, active)


def _start(data, protocols):
    return search(data, protocols, power_shape, (np.arange(.04, .42, .08), np.arange(0., .9, .3)))


def family(name, data, protocols):
    """Returns the shape and its parameter count for one mechanism."""
    if name == 'elastic':
        w = zero(); w[E_COL] = 1.; return w, 0
    if name == 'Hunt-Crossley':
        return free_on(data, protocols, [E_COL, HC]), 1
    if name == 'Kelvin-Voigt':
        return free_on(data, protocols, [E_COL, KV]), 1
    if name == 'Kelvin-Boltzmann, 1 branch':
        best = None
        for i in range(len(TAU)):
            w = free_on(data, protocols, [E_COL, FIRST+i]); j = objective(data, w, protocols)
            if best is None or j < best[0]:
                best = (j, w)
        return best[1], 2
    if name == 'weak power law':
        return search(data, protocols, power_shape, (np.arange(.04, .42, .04), np.arange(0., .9, .15))), 2
    if name == 'weak power law + viscous':
        return search(data, protocols, power_shape, (np.arange(.04, .42, .04), np.arange(0., .9, .15), (0., .0005, .002, .005))), 3
    if name == 'weak power law + Hunt-Crossley':
        def make(b, f, lam):
            w = power_shape(b, f); w[HC] = lam
            return w
        return search(data, protocols, make, (np.arange(.04, .42, .04), np.arange(0., .9, .15), (0., 1e-6, 1e-5, 1e-4))), 3
    if name == 'free spectrum':
        return free_on(data, protocols, np.arange(0, KV), _start(data, protocols)), len(TAU)
    if name == 'free spectrum + viscous':
        return free_on(data, protocols, np.arange(0, KV+1), _start(data, protocols)), len(TAU)+1
    if name == 'free spectrum + Hunt-Crossley':
        return free_on(data, protocols, list(range(0, KV))+[HC], _start(data, protocols)), len(TAU)+1
    if name == 'one branch per decade':
        return free_on(data, protocols, [E_COL]+[FIRST+nearest(x) for x in 10.**np.arange(-3, 3)], _start(data, protocols)), 6
    if name == 'plant hc form':
        return free_on(data, protocols, [E_COL]+[FIRST+nearest(x) for x in (.008, .08, .8)]+[HC]), 4
    if name == 'plant power form':
        return free_on(data, protocols, [E_COL]+[FIRST+nearest(x) for x in (.003, .03, .3)]), 3
    raise ValueError(name)


FAMILIES = ('elastic', 'Hunt-Crossley', 'Kelvin-Voigt', 'Kelvin-Boltzmann, 1 branch', 'plant hc form', 'plant power form', 'weak power law', 'weak power law + viscous', 'weak power law + Hunt-Crossley',
            'one branch per decade', 'free spectrum', 'free spectrum + viscous', 'free spectrum + Hunt-Crossley')
HELD = ('Hunt-Crossley', 'Kelvin-Voigt', 'Kelvin-Boltzmann, 1 branch', 'plant hc form', 'plant power form', 'weak power law', 'weak power law + viscous', 'one branch per decade', 'free spectrum')
ALL = ('rate', 'relax', 'cyclic', 'dma')
RECORDS = None


def one_job(a):
    tissue, name, train = a
    data = Data([r for r in RECORDS if r['tissue'] == tissue])
    w, npar = family(name, data, train)
    return dict(tissue=tissue, family=name, train=list(train), npar=npar, w=w.tolist(), objective=objective(data, w, ALL), score=score(data, w))


def init():
    global RECORDS
    RECORDS = [r for r in prepare() if usable(r)]


def usable(r):
    """Judged on the curves, not on the instrument's flag."""
    if r['kind'] == 'dma':
        return True
    return r['peak_force'] > .3 and r['depth_peak'] > 15. and r['fit_error'] < .08


def main():
    data = prepare()
    print('records', len(data), 'used', sum(usable(r) for r in data))
    jobs = [(t, n, ALL) for t in ('white', 'grey') for n in FAMILIES]
    for t in ('white', 'grey'):
        for n in HELD:
            for hold in ALL:
                jobs.append((t, n, tuple(p for p in ALL if p != hold)))
    with ProcessPoolExecutor(max_workers=14, initializer=init) as ex:
        rows = list(ex.map(one_job, jobs, chunksize=1))
    (HERE/'results'/'neumann_fits.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')
    head = ('rate_nrmse_pct', 'relax_nrmse_pct', 'cyclic_nrmse_pct', 'dma_phase_rmse_deg', 'dma_phase_error_05Hz_deg', 'dma_phase_error_20Hz_deg', 'dma_gain_rmse_pct')
    print('\nall four protocols in the fit: tissue, family, parameters, J, then '+', '.join(head))
    for r in rows:
        if tuple(r['train']) == ALL:
            print(f"{r['tissue'][:5]:5s} {r['family']:32s} {r['npar']:2d} J {1e4*r['objective']:8.1f} "+' '.join(f"{r['score'][k]:7.2f}" for k in head))
    print('\nheld out: the protocol not in the training set is predicted with only the scale fitted')
    for r in rows:
        train = tuple(r['train'])
        if train == ALL:
            continue
        held = [p for p in ALL if p not in train]
        print(f"{r['tissue'][:5]:5s} {r['family']:28s} "+' '.join(f"{p}: {r['score'][p+'_nrmse_pct']:.2f}%" if p != 'dma' else f"dma: phase {r['score']['dma_phase_rmse_deg']:.2f} deg (20 Hz {r['score']['dma_phase_error_20Hz_deg']:+.1f}), gain {r['score']['dma_gain_rmse_pct']:.1f}%" for p in held))


if __name__ == '__main__':
    main()
