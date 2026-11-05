"""Read-only audit of the supplied contact data; no Method, network or plant edits.

Run from any directory: python -X utf8 experiments/rossi_failure_study/real_data_audit.py
Requires numpy, scipy, pandas, numba, pyarrow OR fastparquet, pymupdf. XLSX uses OOXML (stdlib).
Only real_data_audit.json is written in the repository. ZIP extraction, figures
and PDF renders use a TemporaryDirectory outside the repository.
All distributions are descriptive, not population confidence intervals.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
import zipfile
import xml.etree.ElementTree as ET

os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')
_NUMBA_TMP = tempfile.TemporaryDirectory(prefix='contact_audit_numba_')
os.environ['NUMBA_CACHE_DIR'] = _NUMBA_TMP.name
import sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares
from scipy.signal import find_peaks
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = ROOT/'literature/manual_inbox/contact_real_data'
import physics_plant3 as P
import scenes3 as S


def stats(x):
    x = np.asarray(x, float).ravel()
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(n=0)
    return dict(n=len(x), min=float(x.min()), p05=float(np.quantile(x, .05)),
                p25=float(np.quantile(x, .25)), median=float(np.median(x)),
                p75=float(np.quantile(x, .75)), p95=float(np.quantile(x, .95)),
                max=float(x.max()), mean=float(x.mean()), sd=float(x.std(ddof=1)) if len(x)>1 else 0.)


def coeff(row):
    return np.array([row[f'p{i}'] for i in range(1, int(row.fitModel[4:])+2)], float)


def at_force(row, force):
    p = coeff(row); end = float(row.fitMaxDispl)
    grid = np.linspace(0., end, 2001); vals = np.polyval(p, grid)-force
    crossings = np.flatnonzero((vals[:-1] <= 0.) & (vals[1:] >= 0.))
    if len(crossings) != 1:
        return None
    i = crossings[0]; d = brentq(lambda x: np.polyval(p, x)-force, grid[i], grid[i+1])
    k = float(np.polyval(np.polyder(p), d))
    return dict(depth_mm=d, k_N_mm=k, local_exponent=d*k/force)


def remus():
    df = pd.read_csv(RAW/'static/Data Sheet 1.CSV')
    rows = []
    for _, r in df.iterrows():
        p = coeff(r); d = np.linspace(0., min(5., r.fitMaxDispl), 501)
        rows.append(dict(region=r.position, subject=r.id, fit_max_mm=float(r.fitMaxDispl),
                         rmse_N=float(r.rmse), min_tangent_0_5_N_mm=float(np.polyval(np.polyder(p), d).min()),
                         F5_N=float(np.polyval(p, 5.)),
                         by_force={str(f): at_force(r, f) for f in (4, 8, 12)}))
    out = dict(row_count=len(df), regions={}, within_between={}, active_relaxed={},
               method='F polynomial in N, d in mm; solve first unique rising F crossing inside fitMaxDispl; k=dF/dd, n=d*k/F; no extrapolation')
    for region in df.position.unique():
        rr = [r for r in rows if r['region']==region and r['subject'].startswith('P')]
        out['regions'][region] = dict(subjects=len(rr), fit_rmse_N=stats([r['rmse_N'] for r in rr]),
            low_depth_negative_tangent_count=sum(r['min_tangent_0_5_N_mm']<0 for r in rr),
            at_force={str(f): {key: stats([r['by_force'][str(f)][key] for r in rr if r['by_force'][str(f)] is not None])
                for key in ('depth_mm', 'k_N_mm', 'local_exponent')} for f in (4,8,12)},
            range_fits={r['subject']:r['by_force'] for r in rows if r['region']==region and not r['subject'].startswith('P')})
    for state in ('R','A'):
        for f in (4,8,12):
            valid_regions=[state+str(i) for i in range(1,6) if next(r for r in rows if r['region']==state+str(i) and r['subject']=='rangeMean')['by_force'][str(f)] is not None]
            groups = [[r['by_force'][str(f)]['k_N_mm'] for r in rows if r['region']==region
                       and r['subject'].startswith('P') and r['by_force'][str(f)] is not None] for region in valid_regions]
            means = [next(r['by_force'][str(f)]['k_N_mm'] for r in rows if r['region']==region and r['subject']=='rangeMean') for region in valid_regions]
            # Equal weighting of regions; population moments within this finite empirical mixture.
            mu = np.mean([np.mean(g) for g in groups])
            var = np.mean([np.mean((np.array(g)-mu)**2) for g in groups])
            vwithin = np.mean([np.var(g) for g in groups]); vmeans = np.var(means)
            out['within_between'][state+str(f)] = dict(valid_regions=valid_regions,mean_curve_k=stats(means), individual_k=stats(np.concatenate(groups)),
                equal_region_mixture_sd=math.sqrt(var), equal_region_within_sd=math.sqrt(vwithin),
                mean_curve_sd=math.sqrt(vmeans), sd_ratio_individual_over_mean=math.sqrt(var/max(vmeans,1e-20)),
                lost_variance_fraction=1.-vmeans/var)
            if len(means)<2:
                out['within_between'][state+str(f)]['sd_ratio_individual_over_mean']=None
                out['within_between'][state+str(f)]['lost_variance_fraction']=None
    for i in range(1,6):
        r = df[(df.position=='R'+str(i)) & (df.id=='rangeMean')].iloc[0]
        a = df[(df.position=='A'+str(i)) & (df.id=='rangeMean')].iloc[0]
        pr, pa = coeff(r), coeff(a)
        out['active_relaxed'][str(i)] = dict(
            same_depth_force_ratio={str(d):float(np.polyval(pa,d)/np.polyval(pr,d)) for d in (5,10,15)},
            same_force_tangent_ratio={str(f):(at_force(a,f)['k_N_mm']/at_force(r,f)['k_N_mm'] if at_force(a,f) and at_force(r,f) else None) for f in (4,8,12)})
    out['unresolved_force_roots'] = sum(v is None for r in rows if r['subject'].startswith('P') for v in r['by_force'].values())
    out['individual_curve_summaries'] = [r for r in rows if r['subject'].startswith('P')]
    return out


def dahmani():
    a = pd.read_csv(RAW/'static/Acquisition informations.tab', sep='\t', skiprows=1).to_numpy(float)
    a = a[np.isfinite(a[:,1])]
    d = a[:,1]; v = a[:,2::2]; f = a[:,3::2]; df = f-f[0]; rows=[]
    for j in range(6):
        valid=np.isfinite(df[:,j]); dj=d[valid]; fj=df[valid,j]
        # Offset subtraction is a conditional diagnostic, not a recovered contact zero.
        p = np.polyfit(dj, fj, 2); fit = np.polyval(p,dj)
        m = (d>=2.) & (df[:,j]>.5)
        power = np.polyfit(np.log(d[m]),np.log(df[m,j]),1)
        shifted=least_squares(lambda p:p[0]*np.maximum(dj-p[2],0.)**p[1]-fj,
            [.2,2.,3.],bounds=([0.,.2,-3.],[100.,6.,6.]))
        local={str(x):float(np.polyfit(dj[abs(dj-x)<=.61],fj[abs(dj-x)<=.61],1)[0]) for x in (2,4,6)}
        rows.append(dict(acquisition=j+1,valid_levels=int(valid.sum()),max_depth_mm=float(dj[-1]),
            offset_N=float(f[0,j]),delta_F_at_last_N=float(fj[-1]),delta_F7_N=float(df[-1,j]) if valid[-1] else None,
            quadratic_coeff_N_mm=p.tolist(), quadratic_residual_rms_N=float(np.sqrt(np.mean((fj-fit)**2))),
            tangent_N_mm={str(x):float(np.polyval(np.polyder(p),x)) for x in (2,4,6)},
            local_linear_tangent_N_mm=local,
            log_power_2_7mm=float(power[0]), power_points=int(m.sum()),
            shifted_power=dict(a=float(shifted.x[0]),exponent=float(shifted.x[1]),depth_origin_mm=float(shifted.x[2]),
                residual_rms_N=float(np.sqrt(np.mean(shifted.fun**2)))),
            apparent_scatter_0_3mm_N=float(np.std(fj[dj<=3.],ddof=1))))
    unique_v=np.unique(v[np.isfinite(v)]); unique_f=np.unique(f[np.isfinite(f)])
    out=dict(levels=len(d), acquisitions=6, raw_force_N=stats(f), rows=rows,
        conditional_increment_N_at_2_4_6p4_7={str(x):stats(df[np.argmin(abs(d-x))]) for x in (2,4,6.4,7)},
        voltage_grid_min_step_V=float(np.diff(unique_v).min()), force_level_step_N=stats(np.diff(unique_f)),
        raw_force_level_count=len(unique_f),within_depth_across_location_sd_N=stats(np.nanstd(f,axis=1,ddof=1)),
        increment_across_location_sd_N=stats(np.nanstd(df,axis=1,ddof=1)),
        caveat='six acquisitions at different locations, not six technical repeats; residual contains shape error and quantization, not identified sensor noise; offset cause unidentifiable')
    # Quadratic V->F calibration can itself create force curvature.
    valid=np.isfinite(v)&np.isfinite(f);q=np.polyfit(v[valid],f[valid],2)
    out['voltage_force_quadratic']=dict(coeff=q.tolist(),max_residual_N=float(np.max(abs(np.polyval(q,v[valid])-f[valid]))))
    return out


def xlsx_tables(path):
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        strings=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        names=[e.attrib['name'] for e in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets',ns)]
        tables={}
        for i,name in enumerate(names,1):
            rr=[]
            for row in ET.fromstring(z.read(f'xl/worksheets/sheet{i}.xml')).find('m:sheetData',ns):
                dd={}
                for c in row:
                    v=c.find('m:v',ns)
                    if v is not None:
                        dd[re.sub(r'\d','',c.attrib['r'])]=strings[int(v.text)] if c.attrib.get('t')=='s' else float(v.text)
                rr.append(dd)
            tables[name]=rr
        return tables


def callera():
    tabs=xlsx_tables(RAW/'dynamic/repository_test_indentation.xlsx')
    dma=[r for r in tabs['DMA'][1:] if isinstance(r.get('B'),float)]
    creep=tabs['creep'][1:]; areas=tabs['hysteresis area'][1:]
    pause=tabs['pause'][2:]; conditioning={}
    for label,col,testcol in [('pause','C','B'),('no_pause','G','F')]:
        triples=[]
        for k in range(0,len(pause)-2,3):
            g=pause[k:k+3]
            if [r.get(testcol) for r in g]==[1.,2.,3.] and all(isinstance(r.get(col),float) for r in g):
                triples.append(g[2][col]/g[0][col])
        conditioning[label+'_third_over_first']=stats(triples)
    ep=np.array([r['B'] for r in dma]); epp=np.array([r['C'] for r in dma]); tan=np.array([r['D'] for r in dma])
    return dict(sheet_rows={k:len(v) for k,v in tabs.items()},
        header_units={k:(v[:2] if k=='pause' else v[:1]) for k,v in tabs.items()},
        DMA_Eprime_kPa=stats(ep),DMA_Edoubleprime_kPa=stats(epp),
        DMA_reported_tandelta=stats(tan),DMA_ratio_Edouble_over_Eprime=stats(epp/ep),
        DMA_reported_phase_deg=stats(np.degrees(np.arctan(tan))),
        DMA_quotient_phase_deg=stats(np.degrees(np.arctan(epp/ep))),
        DMA_max_tan_discrepancy=float(np.max(abs(tan-epp/ep))),
        creep_status_counts={label:sum(r.get('D')==label for r in creep) for label in sorted(set(r.get('D','') for r in creep))},
        creep_plateau_s=stats([r['C'] for r in creep if r.get('D')=='Plateau' and 'C' in r]),
        creep_force_column_as_labelled_N=stats([r['B'] for r in creep if 'B' in r]),
        hysteresis_area_microN_microm=stats([r['C'] for r in areas if 'C' in r]),
        conditioning=conditioning,
        limitations='summary workbook only; E_Hz is a column label, not frequency; DMA frequency/amplitude and indentation ramp/time-series absent; creep Force [N] label not independently verified; no energy denominator')


def runs(mask):
    edges=np.diff(np.r_[False,mask,False].astype(int))
    return list(zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1)-1))


def relaxation_metrics(t,y):
    """Descriptive, positive exponential fits on log-time interpolation, not mechanism ID."""
    t=np.asarray(t)-t[0]; y=np.asarray(y)/y[0]
    end=min(float(t[-1]),290.)
    grid=np.geomspace(1.,end,100); obs=np.interp(grid,t,y)
    one=least_squares(lambda p:p[0]+p[1]*np.exp(-grid/p[2])-obs,
                      [.4,.6,20.],bounds=([0,0,.001],[2,2,10000]))
    best=None
    for guess in ([.3,.4,.3,2.,100.],[.3,.3,.4,10.,1000.]):
        two=least_squares(lambda p:p[0]+p[1]*np.exp(-grid/p[3])+p[2]*np.exp(-grid/p[4])-obs,
                          guess,bounds=([0,0,0,.001,.001],[2,2,2,10000,10000]),max_nfev=2000)
        if best is None or np.mean(two.fun**2)<np.mean(best.fun**2):best=two
    knots=[1.,3.,10.,30.,100.,end]
    values=np.interp(knots,t,y)
    return dict(duration_s=float(t[-1]),normalized_at_s=dict(zip(map(str,knots),values.tolist())),
        end_over_start=float(np.median(y[t>=t[-1]-5.])),
        log_time_slopes_per_decade={f'{a:g}-{b:g}':float((np.interp(b,t,y)-np.interp(a,t,y))/math.log10(b/a)) for a,b in [(1,10),(10,100),(100,end)]},
        single_exp=dict(offset=float(one.x[0]),amplitude=float(one.x[1]),tau_s=float(one.x[2]),rms=float(np.sqrt(np.mean(one.fun**2)))),
        two_exp=dict(parameters_c_a_b_tau1_tau2=best.x.tolist(),rms=float(np.sqrt(np.mean(best.fun**2)))))


def cycle_work(d,f):
    inc=np.diff(d)*.5*(f[:-1]+f[1:]); win=float(np.maximum(inc,0).sum()); loss=float(inc.sum())
    return dict(input_work=win,net_work=loss,loss_over_input=loss/win if win>0 else None)


def brain():
    out={}; protocols=[]
    for region,height in [('corona_radiata',5.),('visual_cortex',3.4)]:
        a=np.loadtxt(RAW/f'dynamic/cyclic_relaxation_{region}.csv',delimiter=',')
        t,f,d=a.T; amp=float(np.max(abs(d))); holds=[]
        # A 0.5% displacement tolerance detects plateaus; retain only >200 s.
        for start,end in runs(abs(d)>.995*amp):
            if t[end]-t[start]>200.:
                # Avoid the final fraction of the loading ramp; report the 0.5 s exclusion.
                idx=np.flatnonzero((t>=t[start]+.5)&(t<=min(t[end],t[start]+.5+290.)))
                yy=f[idx]*np.sign(d[idx[0]])
                holds.append(dict(sign='compression' if d[start]<0 else 'tension',start_s=float(t[idx[0]]),
                    displacement_mm=float(np.median(d[idx])),displacement_span_mm=float(np.ptp(d[idx])),
                    metrics=relaxation_metrics(t[idx],yy)))
        neg=find_peaks(-d,prominence=.5*amp,distance=70)[0]
        pos=find_peaks(d,prominence=.5*amp,distance=70)[0]
        neg=[i for i in neg if t[i]<holds[0]['start_s']-10.]
        pos=[i for i in pos if t[i]<holds[0]['start_s']-10.]
        # Merge duplicate positive turns in the data-reduced plateau (time separation <10 s).
        pos=[i for j,i in enumerate(pos) if j==0 or t[i]-t[pos[j-1]]>10.]
        cyc=[]; start=0
        for ni,pi in zip(neg[:3],pos[:3]):
            end_candidates=np.flatnonzero((np.arange(len(d))>pi)&(d<=.002*amp))
            end=int(end_candidates[0]); tt=t[start:end+1]; dd=d[start:end+1]; ff=f[start:end+1]
            # Near-zero displacement crossing; endpoint residual is <0.2% amplitude.
            work=cycle_work(dd,ff)
            cyc.append(dict(duration_s=float(tt[-1]-tt[0]),amplitude_mm=float(max(abs(dd))),
                compression_peak_force_N=float(-ff.min()),tension_peak_force_N=float(ff.max()),
                work_N_mm=work))
            start=end
        # First unilateral compression excursion: used for matched shape-only predictions.
        ni=neg[0]; end=int(np.flatnonzero((np.arange(len(d))>ni)&(d>=0.))[0])
        protocols.append((region,'first_compression',t[:end+1],np.maximum(-d[:end+1]/amp,0.),np.maximum(-f[:end+1],0.)))
        hh=holds[0]; hidx=int(np.searchsorted(t,hh['start_s'])); loadstart=max(0,int(np.flatnonzero((t<t[hidx])&(abs(d)<.01*amp))[-1]))
        protocols.append((region,'compression_hold',t[loadstart:np.searchsorted(t,hh['start_s']+290.)],
                          np.maximum(-d[loadstart:np.searchsorted(t,hh['start_s']+290.)]/amp,0.),None))
        out[region]=dict(rows=len(t),duration_s=float(t[-1]),
            interpreted_columns='time s, force N, displacement mm (CSV has no header; identification from paper geometry and loading velocity)',
            height_mm=height,peak_strain=float(amp/height),time_step_s=stats(np.diff(t)),
            force_N=stats(f),cycles=cyc,holds=holds,
            first_compression_work_N_mm=cycle_work(-d[:end+1],-f[:end+1]),
            preprocessing='published moving average + RDP point reduction; small minimum dt is not usable acquisition bandwidth')
    return out,protocols


def rheometer(tmp):
    path=RAW/'dynamic/Rheometer data.zip'
    with zipfile.ZipFile(path) as z:
        # Check every resolved destination before extracting untrusted archive members.
        dest=tmp/'rheometer'; dest.mkdir()
        for member in z.infolist():
            target=(dest/member.filename).resolve()
            if not target.is_relative_to(dest.resolve()):raise ValueError('ZIP path escape')
        z.extractall(dest)
    entries=[]; geometries=[]; duration=[]; raw_counts=[]
    for g in sorted(dest.rglob('geometry.yaml')):
        text=g.read_text(); h=float(re.search(r'height:\s*([\deE.+-]+)',text)[1]); geometries.append(h*1000.)
        for sign in ('neg','pos'):
            group=[]
            for k in (1,2,3):
                a=pd.read_csv(g.parent/f'comp_ten_cyc_l1_c{k}_{sign}.csv',comment='#')
                x=a.displacement.to_numpy(float); f=a.force.to_numpy(float)
                peak=float(np.max(abs(f))); amp=float(np.max(abs(x))/h)
                group.append(dict(cycle=k,processed_peak_force_N=peak,strain=amp,
                    displacement_return_fraction=float(abs(x[-1]-x[0])/max(np.ptp(x),1e-20)),
                    time_is_row_index=bool(np.array_equal(a.time,np.arange(len(a))))))
                duration.append(float(a.time.iloc[-1]-a.time.iloc[0])); raw_counts.append(len(a))
            entries.append(dict(specimen=g.parent.name,sign=sign,cycles=group,
                third_over_first_peak=group[2]['processed_peak_force_N']/group[0]['processed_peak_force_N']))
    return dict(zip_md5=hashlib.md5(path.read_bytes()).hexdigest(),specimens=len(geometries),
        cord_groups=len(list(dest.glob('Rheometer data/Spinal cord *'))),height_mm=stats(geometries),
        per_curve_points=stats(raw_counts),time_column_span=stats(duration),
        hysteresis_loss_ratio=None,closed_cycle_count=sum(c['displacement_return_fraction']<.05 for e in entries for c in e['cycles']),
        time_index_count=sum(c['time_is_row_index'] for e in entries for c in e['cycles']),
        third_over_first_peak={s:stats([e['third_over_first_peak'] for e in entries if e['sign']==s]) for s in ('neg','pos')},
        strain=stats([c['strain'] for e in entries for c in e['cycles']]),entries=entries,
        provenance='MD5 matches Zenodo 17152951, porcine spinal cord multimodal pipeline; not 14536174 indentation archive',
        units='geometry in m (radius .004); displacement m and force N consistent with SI preprocessor; time is integer row index 0..n-1 after reduction, not seconds; no physical rate or Hz inferred',
        limitations='hyperelastic preprocessing averages loading/unloading; individual c1/c2/c3 files also have one-way displacement branches, not closed loops; neither those nor mean_neg/pos can measure hysteresis')


def simulate_material(t,d,alpha,tau,lam):
    q=np.zeros(3); result=np.zeros(len(t)); result[0]=d[0]
    for i in range(1,len(t)):
        dt=t[i]-t[i-1]
        if dt<=0:raise ValueError('nonmonotonic protocol')
        e=np.exp(-dt/tau); delta=d[i]-d[i-1]
        q=q*e+alpha*tau/dt*(1.-e)*delta
        result[i]=max((d[i]+q.sum())*(1.+lam*delta/dt),0.)
    return result


def dynamics(protocols):
    freqs=np.array([.01,.1,1.,5.,10.,20.,20.8,60.]); rows=[]
    for kind in ('hc','power','kb'):
        for target in (7.,11.,17.,23.):
            for seed in range(3):
                alpha,tau,phi=S._dynamic_response(kind,target,np.random.default_rng(seed))
                lam=math.tan(math.radians(phi))/(2.*math.pi*20.8)
                z=2j*np.pi*freqs[:,None]*tau
                H=1.+np.sum(alpha*z/(1.+z),axis=1)+2j*np.pi*freqs*lam
                relaxation=(1.+np.exp(-np.array([0.,.01,.1,1.,10.,100.,300.])[:,None]/tau).dot(alpha))/(1.+sum(alpha))
                matches=[]
                for region,label,t,d,f in protocols:
                    y=simulate_material(t,d,alpha,tau,lam)
                    if label=='first_compression':
                        matches.append(dict(region=region,protocol=label,model_loss_ratio=cycle_work(d,y)['loss_over_input'],
                            measured_loss_ratio=cycle_work(d,f)['loss_over_input']))
                    else:
                        # Sustained compression begins when normalized d reaches .995.
                        ii=int(np.flatnonzero(d>.995)[0]); jj=int(np.searchsorted(t,t[ii]+1.))
                        matches.append(dict(region=region,protocol=label,hold_end_over_1s=float(y[-1]/y[jj])))
                rows.append(dict(kind=kind,target_phase_deg=target,seed=seed,alpha=alpha.tolist(),tau_s=tau.tolist(),
                    phi_hc_deg=phi,lambda_normalized=lam,phase_deg=np.degrees(np.angle(H)).tolist(),
                    modulus_over_equilibrium=abs(H).tolist(),step_normalized=relaxation.tolist(),protocols=matches))
    return dict(frequencies_Hz=freqs.tolist(),step_times_s=[0.,.01,.1,1.,10.,100.,300.],rows=rows,
        convention='analytic frozen material H=1+sum(alpha*iw*tau/(1+iw*tau))+iw*lambda, F0=k=1; same measured normalized compression displacement/time protocols; equilibrium modulus set to 1, no abdominal parameter transfer; no arm/sensor/support',
        caveat='HC normalized lambda uses F0/k=1 depth unit; large compression protocol is a shape diagnostic, not a fit of glued tensile brain data or finite probe contact')


def scanning():
    a=pd.read_parquet(RAW/'scanning/USForce_episode_000000.parquet',columns=['state','actions','timestamp','stage','prompt'])
    s=np.stack(a.state); action=np.stack(a.actions); t=a.timestamp.to_numpy(float); dt=np.diff(t)
    f=s[:,6:9]; norm=np.linalg.norm(f,axis=1); contact=norm>0
    # Hypothesis from coordinate magnitudes: positions m, angles rad.
    # Card explicitly labels roll/pitch/yaw, so Euler is primary; UR rotvec is an alternative.
    v=np.diff(s[:,:3],axis=0)/dt[:,None]*1000.; speed=np.linalg.norm(v,axis=1)
    acc=np.diff(v,axis=0)/(.5*(dt[1:]+dt[:-1]))[:,None]/1000.
    rv=Rotation.from_rotvec(s[:,3:6]); eu=Rotation.from_euler('xyz',s[:,3:6])
    def rotstat(rot):
        omega=(rot[:-1].inv()*rot[1:]).magnitude()/dt
        axes=rot.apply(np.tile([0.,0.,1.],(len(s),1)))
        axial=np.sum(v*.5*(axes[:-1]+axes[1:]),axis=1)
        projected=np.sum(f*axes,axis=1)
        return dict(angular_rate_deg_s=stats(np.degrees(omega)),axial_velocity_mm_s=stats(axial),
            base_frame_force_projected_tool_z_assumption=stats(projected[contact]),
            orientation_change_deg=float(np.degrees((rot[0].inv()*rot[-1]).magnitude())))
    first=int(np.flatnonzero(contact)[0]); intervals=contact[:-1]&contact[1:]
    st=a.stage.to_numpy(); stable=st==1;stage_intervals=stable[:-1]&stable[1:]
    changes=np.r_[True,np.any(s[1:,6:]!=s[:-1,6:],axis=1)]
    runstarts=np.flatnonzero(changes); runlengths=np.diff(np.r_[runstarts,len(s)])
    # Actions describe changes between successive recorded states; check scales/rotation increments.
    return dict(frames=len(s),duration_s=float(t[-1]-t[0]),dt_s=stats(dt),
        state_column_ranges=[stats(s[:,i]) for i in range(12)],
        initial_zero_force_frames=first,first_nonzero_time_s=float(t[first]),
        all_six_wrench_initially_exact_zero=bool(np.all(s[:first,6:]==0.)),
        stage_counts={str(k):int((a.stage==k).sum()) for k in a.stage.unique()},
        first_nonzero_stage=int(a.stage.iloc[first]),prompt=str(a.prompt.iloc[0]),
        translation_action_vs_difference_max_error=float(np.max(abs(action[:-1,:3]-np.diff(s[:,:3],axis=0)))),
        force_magnitude_if_N=stats(norm[contact]),force_magnitude_all_frames_if_N=stats(norm),
        stage1_force_magnitude_if_N=stats(norm[stable]),
        wrench_update_count=int(changes.sum()),identical_wrench_run_frames=stats(runlengths),
        fz_if_tool_frame_N=stats(f[contact,2]),
        speed_if_mm_s=stats(speed),nonzero_wrench_speed_if_mm_s=stats(speed[intervals]),
        stage1_speed_if_mm_s=stats(speed[stage_intervals]),
        acceleration_if_m_s2=stats(np.linalg.norm(acc,axis=1)),
        nonzero_wrench_acceleration_if_m_s2=stats(np.linalg.norm(acc,axis=1)[intervals[:-1]&intervals[1:]]),
        path_length_if_mm=float(np.linalg.norm(np.diff(s[:,:3],axis=0),axis=1).sum()*1000.),
        rotation_vector_hypothesis=rotstat(rv),euler_xyz_hypothesis=rotstat(eu),
        inference='meters and radians plausible from 0.48..0.58 coordinate and UR-sized workspace, not specified in data card; Euler follows roll/pitch/yaw labels, rotvec is alternative. Global pose frame, force frame, calibration, TCP/probe axis and force sign unresolved. Nonzero wrench is not confirmed contact. Magnitude is wrench-force magnitude, not established tissue normal force. Initial exact zero block may be unavailable/uninitialized wrench or noncontact, not noise proof.',
        bandwidth='30 FPS, Nyquist 15 Hz; no 20.8 Hz phase, intrinsic noise, actuator lag, indentation or body motion identification')


def consistency():
    out={}; rows=[]
    # No scene simulation: static contact calls at a fixed coordinate only.
    for law in ('foundation','halfspace'):
        for region in ('R1','R3','A1','A5'):
            sc=S.cell_scene3('rest','relaxed',T=.1,overrides=dict(contact_law=law,region=region,
                probe='convex',slide='long',field_sigma=0.,surface_amp=0.,pose_error=0.,force=8.))
            tab=sc['curve_table'];g=P._geometry(sc)
            for d in (.002,.008,.015,.025):
                f,k,*_=P.contact(sc,d,0.)
                eps=1e-7; fd=(P.contact(sc,d+eps,0.)[0]-P.contact(sc,d-eps,0.)[0])/(2*eps)
                rows.append(dict(law=law,region=region,depth_mm=d*1000.,reported_k_N_mm=k/1000.,
                    finite_difference_k_N_mm=fd/1000.,relative_k_error=k/fd-1.))
            # Uniform full footprint on a thin cover: designed rib cap versus actual mapped slope.
            hard=dict(sc,full_cover=True,cover=.007,k_rib=10000.)
            ratios=[]
            for d in np.linspace(.007,.05,50):
                ep=1e-7;fd=(P.contact(hard,d+ep,0.)[0]-P.contact(hard,d-ep,0.)[0])/(2*ep)
                ratios.append(fd/hard['k_rib'])
            rows[-1]['rib_cap_mapped_tangent_over_krib_range']=[float(min(ratios)),float(max(ratios))]
            rows[-1]['table_last_map_slope']=float(tab[2][-1])
    out['static_mapping']=rows
    # Scene total target is material phase, not series or measured phase.
    phases=[]
    for kind in ('hc','power','kb'):
        for target in (7.,11.,17.,23.):
            sc=S.cell_scene3('rest','relaxed',T=.1,overrides=dict(dynamics=kind,phase=target,
                region='R1',field_sigma=0.,surface_amp=0.,pose_error=0.,force=8.))
            ref=dict(sc,q_field=np.full_like(sc['q_field'],sc['q_base']),E_field=np.full_like(sc['E_field'],sc['E_base']))
            d=P.solve_depth(ref,8.,float(sc['xcmd'][0])); k=P.contact(ref,d)[1]; w=2*np.pi*20.8
            H=k*(1+sum(sc['alpha']*1j*w*sc['tau']/(1+1j*w*sc['tau'])))+1j*w*P.calibration(sc)*8.
            sensor_pole=math.exp(-2*math.pi*30*.001)
            sensor_phase=float(np.degrees(np.angle((1.-sensor_pole)/(1.-sensor_pole*np.exp(-1j*w*.001)))))
            phases.append(dict(kind=kind,target=target,material_deg=float(np.degrees(np.angle(H))),
                series_deg=float(np.degrees(np.angle(P.complex_gain(sc,k,8.)))),
                sensor_IIR_phase_deg=sensor_phase,
                series_plus_sensor_IIR_deg=float(np.degrees(np.angle(P.complex_gain(sc,k,8.))))+sensor_phase))
    out['phase_labels']=phases
    sensor_pole=math.exp(-2*math.pi*30*.001)
    out['noise']=dict(input_sigma_N=.08,output_sigma_N_30Hz=.08*math.sqrt((1-sensor_pole)/(1+sensor_pole)),
                      bias_random_walk_sd_10s_N=.0002*math.sqrt(10000.),
                      readme='sensor sigma is injected before IIR; random walk is separate')
    # Actual aggregate drive velocity in representative combination, no dynamics/Method.
    drive=[]
    for i in range(4):
        sc=S.cell_scene3('combo','relaxed',i=i,T=10.,overrides=dict(no_curve=True))
        drive.append(float(np.max(abs(np.diff(sc['u'])))/.001))
    out['combo_drive_peak_m_s']=drive
    transfers=[]
    means=pd.read_csv(RAW/'static/Data Sheet 1.CSV')
    for region,poly in S.REMUS.items():
        maxd=float(means[(means.position==region)&(means.id=='rangeMean')].fitMaxDispl.iloc[0])
        for probe,geo in S.PROBES.items():
            for gamma in (.5,1.):
                depths=np.geomspace(1.e-5,.06,1000);F,k=P.probe_curve(np.array(poly),geo,gamma,depths)
                aa=np.minimum(np.sqrt(2.*geo['radius']*depths),geo['length']/2.) if math.isfinite(geo['radius']) else np.full_like(depths,geo['length']/2.)
                area=geo['width']*2.*aa;dh=np.geomspace(1.e-5,.05,4000)
                equiv=np.maximum(np.interp(depths/np.sqrt(area),dh/np.sqrt(P._hemisphere_area(dh)),dh),.005)*1000.
                transfers.append(dict(region=region,probe=probe,gamma=gamma,reference_fit_max_mm=maxd,
                    minimum_tangent_N_mm=float(k.min()/1000.),negative_force_increment_count=int(np.sum(np.diff(F)<0)),
                    first_negative_tangent_depth_mm=float(depths[np.flatnonzero(k<0)[0]]*1000.) if np.any(k<0) else None,
                    by_force={str(f):dict(probe_depth_mm=float(np.interp(f,F,depths)*1000.),
                        equivalent_reference_depth_mm=float(np.interp(f,F,equiv)),
                        uses_clamped_toe=bool(np.interp(f,F,equiv)<=5.000001),
                        extrapolates=bool(np.interp(f,F,equiv)>maxd)) for f in (4,8,12)} if np.all(np.diff(F)>=0) else None))
    out['transfer_domains']=transfers
    provenance=[]
    for region,poly in S.REMUS.items():
        r=means[(means.position==region)&(means.id=='rangeMean')].iloc[0];p=coeff(r)
        provenance.append(dict(region=region,code_degree=len(poly)-1,S1_degree=len(p)-1,
            depth_mm=[5,10,20],code_force_N=np.polyval(poly,[5,10,20]).tolist(),
            S1_force_N=np.polyval(p,[5,10,20]).tolist(),
            code_force_at_fitmax_N=float(np.polyval(poly,r.fitMaxDispl)),S1_force_at_fitmax_N=float(np.polyval(p,r.fitMaxDispl))))
    out['code_vs_S1_coefficients']=provenance
    qscale=[]
    for scale in (.5,1.,2.):
        sc=S.cell_scene3('rest','relaxed',T=.1,overrides=dict(q_scale=scale,region='R1',
            field_sigma=0.,surface_amp=0.,pose_error=0.,probe='linear',slide='long',contact_law='foundation'))
        qscale.append(dict(q_scale=scale,F_at_15mm_N=float(P.contact(sc,.015,0.)[0])))
    out['q_scale_cancellation']=qscale
    a2=S.cell_scene3('rest','active',T=.1,overrides=dict(region='A2',probe='linear',
        contact_law='foundation',field_sigma=0.,surface_amp=0.,pose_error=0.))
    ga2=P._geometry(a2)
    out['A2_last_slope_cap_hazard']=dict(last_map_slope=float(a2['curve_table'][2][-1]),
        foundation_cap_N_m3=float(ga2[13]),halfspace_cap_Pa=float(ga2[14]),
        interpretation='negative last map slope makes internal cover caps negative; consequence of inconsistent A2 source and unchecked extrapolation, not a measured bone property')
    return out


def beber(tmp):
    import pymupdf as fitz
    doc=fitz.open(ROOT/'literature/papers/beber.pdf')
    # PyMuPDF and Cairo emit different SVG coordinates. Use PDF drawing paths directly.
    candidates=[]
    for drawing in doc[3].get_drawings():
        col=drawing.get('color')
        pts=[]
        for item in drawing['items']:
            if item[0]=='l':pts.extend([(item[1].x,item[1].y),(item[2].x,item[2].y)])
        if col is not None and len(pts)>1000:
            candidates.append((drawing,np.asarray(pts)))
    # Main blue curve: stroke near (0,.447,.741), horizontal axis span >150 PDF points.
    blue=[(dr,a) for dr,a in candidates if dr['color'][2]>.5 and dr['color'][0]<.2 and np.ptp(a[:,0])>150.]
    if not blue:raise ValueError('Beber vector curve not found')
    dr,xy=min(blue,key=lambda v:np.ptp(v[1][:,0]));
    # Main PDF 0/15 ticks in points; ignore the outer border and rotated inset.
    dd=(xy[:,0]-334.1180725)/(553.118042-334.1180725)*15.
    ff=(136.3317108-xy[:,1])/(136.3317108-53.49012375)*15.
    def smooth(x):return np.convolve(np.pad(x,20,mode='edge'),np.ones(41)/41,mode='valid')
    dd,ff=smooth(dd),smooth(ff); top=int(np.argmax(dd))
    le=int(np.flatnonzero(dd>dd.max()-.3)[0]); ue=top
    grid=np.arange(1.,14.1,1.)
    order=np.argsort(dd[:le]); order2=np.argsort(dd[ue:])
    fl=np.interp(grid,dd[:le][order],ff[:le][order]);fu=np.interp(grid,dd[ue:][order2],ff[ue:][order2])
    loop=float(np.trapz(fl-fu,grid)/np.trapz(fl,grid))
    from check_curve3 import BEBER_D,BEBER_F
    return dict(vector_blue_segment_endpoints=len(xy),
        partial_loop_1_14mm_loss_over_loading_work=loop,vector_peak_depth_mm=float(dd.max()),vector_peak_force_N=float(ff.max()),
        vector_loading_exponent_1_14mm=float(np.polyfit(np.log(grid),np.log(fl),1)[0]),
        vector_loading_tangent_N_mm=stats(np.gradient(fl,grid)),
        tabulated_loading_exponent_1_14mm=float(np.polyfit(np.log(BEBER_D),np.log(BEBER_F),1)[0]),
        tabulated_loading_tangent_N_mm=stats(np.gradient(BEBER_F,BEBER_D)),
        note='Vector re-extracted from local PDF; partial loop 1..14mm, not complete closed-cycle loss. No time base in force-depth figure; hysteresis must not be converted to 20.8Hz phase.')


def main():
    with tempfile.TemporaryDirectory(prefix='real_data_audit_') as temp:
        tmp=Path(temp)
        b,protocols=brain()
        out=dict(scope='analysis only; no Method; no plant or scene edits; no data downloads',
            input_sha256={str(p.relative_to(RAW)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(RAW.rglob('*')) if p.is_file()},
            source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in (HERE/'physics_plant3.py',HERE/'scenes3.py',HERE/'check_curve3.py',
                                     HERE/'beber_fig2.py',ROOT/'literature/papers/beber.pdf')},
            remus=remus(),dahmani=dahmani(),callera=callera(),brain=b,rheometer=rheometer(tmp),
            dynamics=dynamics(protocols),scanning=scanning(),consistency=consistency(),beber=beber(tmp))
        target=HERE/'results/real_data_audit.json'
        target.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print('Wrote',target)
        print('Remus within/between:',{k:(round(v['sd_ratio_individual_over_mean'],3) if v['sd_ratio_individual_over_mean'] is not None else None) for k,v in out['remus']['within_between'].items()})
        print('Brain hold ratios:',{k:[round(h['metrics']['end_over_start'],3) for h in v['holds']] for k,v in b.items()})
        print('USForce nonzero-wrench speed:',out['scanning']['nonzero_wrench_speed_if_mm_s'])


if __name__=='__main__':
    main()
