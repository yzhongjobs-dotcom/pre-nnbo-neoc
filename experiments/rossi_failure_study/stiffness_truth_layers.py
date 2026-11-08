"""Noise-free open-loop layer audit; existing plant equations are unchanged.

python -X utf8 experiments/rossi_failure_study/stiffness_truth_layers.py
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import sys
sys.dont_write_bytecode=True
import copy
import hashlib
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import low_estimate_forensics as L
import corrected_comparison as C

P=L.P3
HERE=Path(__file__).resolve().parent
OUT=HERE/'results/stiffness_truth_layers.json'
SOURCE=HERE/'results/id_layer_root_cause.json'
FREQUENCIES=(5.,10.,15.,20.8,30.,35.)
LABELS=('A_neural_white_i6','A_slow_thin_neural_white_i1','A_slow_thin_viscous_i5')
DURATION=16.
START=8.
TC=.008
DT=.0005


def describe(z):
    z=complex(z)
    return dict(real=z.real,imag=z.imag,gain=abs(z),phase_deg=float(np.angle(z)*180/np.pi))


def phasor(t,v,f):
    s=t-t.mean()
    A=np.c_[np.cos(2*np.pi*f*t),np.sin(2*np.pi*f*t),np.ones(len(t)),s,s*s]
    b=np.linalg.lstsq(A,np.asarray(v,float),rcond=None)[0]
    return complex(b[0],-b[1])


def summary(a):
    a=np.asarray(a,float);a=a[np.isfinite(a)]
    return dict(n=len(a),median=float(np.median(a)),min=float(a.min()),max=float(a.max())) if len(a) else dict(n=0)


def scene(row):
    """Static original spatial field at the recorded t=4 s scan coordinate.

    No tissue parameters are changed, no travel, no time-varying loads. Fully
    relaxed Prony state is initialized exactly, rather than waiting 25 s tails.
    Infinitesimal lateral motion remains part of the original mechanics and
    its excursion is measured, as is the local K variation under the sine.
    """
    L.O.setup()
    sc=L.MC.scene_for(tuple(row['case']['cell']),C.BASE,row['case']['variant'])
    old=np.asarray(row['records']['log'],float)
    j=int(np.argmin(abs(old[:,0]-4.)))
    x=float(old[j,1]);size=int(DURATION*1000)+10
    sc=copy.deepcopy(sc);sc['T']=DURATION
    for key in ('xcmd','u','fext','act','rip'):
        sc[key]=np.full(size,x if key=='xcmd' else 0.)
    sc['memory_start']=np.zeros_like(sc['alpha'])
    sc['start_gap']=0.
    return sc,dict(source_time=float(old[j,0]),fixed_x_command=x,
                   previous_moving_gain_ratio=float(np.exp(row['sampled_gap']['bias'])))


def initial(sc,intervention):
    sc=copy.deepcopy(sc)
    if intervention=='mass0':sc['tool_mass']=sc['tool_mass_hat']=0.
    elif intervention=='kr_large':sc['kr']=1.e12
    elif intervention=='offset0':sc['sensor_delay_ms']=sc['sensor_delay_hat']=0.
    pl=dict(P.PLATFORM,tau_a=sc.get('tau_a_true',P.PLATFORM['tau_a']),log_dt=DT)
    s,buf=P._init(sc,pl);c=P._constants(sc,pl)
    if intervention=='fast1ms':c[35]=.001
    elif intervention=='fast4ms':c[35]=.004
    return sc,pl,s,buf,c


def constitutive(sc,state,c,frequency):
    g=P._geometry(sc);x=state[7];d=state[13];theta=state[17]
    fe,k,cop,length,cover=P._contact(x,d,theta,0.,g,sc['q_field'],sc['E_field'],sc['ribs'],0.)
    # Exactly the logger's derivative convention, separately from its fixed f.
    if g[5]!=0.:
        k=(P._contact(x,d+1e-6,theta,0.,g,sc['q_field'],sc['E_field'],sc['ribs'],0.)[0]-
           P._contact(x,d-1e-6,theta,0.,g,sc['q_field'],sc['E_field'],sc['ribs'],0.)[0])/2e-6
    w=2*np.pi*frequency;alpha,tau,_=P.memory_arrays(sc)
    kc=k*(1+np.sum(alpha*(1j*w*tau)/(1+1j*w*tau)))+1j*w*c[4]*fe
    ct=np.cos(theta);st=np.sin(theta)
    fr=P._friction(1.,0.,length,c,g)-c[17]*0.
    axial=ct+fr*st;lateral=fr*ct-st
    series=ct*(axial/c[3]+c[26]*lateral)
    return kc,kc*axial*ct,kc*axial*ct/(1+kc*series),dict(
        normal_elastic_tangent=float(k),elastic_force=float(fe),series=float(series),
        theta_deg=float(theta*180/np.pi),depth=float(d),normal_force=float(state[14]))


def run(sc0,frequency,intervention='base',phase=0.,amplitude=None,check=False):
    sc,pl,s,buf,c=initial(sc0,intervention)
    start_state=s.copy();kc,ka,kt,op=constitutive(sc,s,c,frequency)
    # Target roughly 0.12 N fundamental at the contact, still small versus 8 N.
    # A single amplitude is held across interventions for a given condition/f.
    w=2*np.pi*frequency
    servo=1/((1+1j*w*(pl['tau_a']-.002))*(1+1j*w*.002))
    if amplitude is None:amplitude=min(2e-4,.12/max(abs(kt*servo),1.))
    ncmd=int(round(DURATION/TC));tt=np.arange(ncmd)*TC
    command=float(s[2])+amplitude*np.sin(w*tt+phase)
    alpha,tau,_=P.memory_arrays(sc);g=P._geometry(sc)
    noise=np.zeros((int(DURATION*1000)+10,6))
    log=np.zeros((ncmd*16+2,30));n=kk=0
    # Diagnostic tick adds xc/acceleration columns; all force equations are
    # identical to P._tick, with its optional spatial freeze flag kept zero.
    tick=L.diagnostic_tick()
    for cmd in command:
        seen=s[5];n0=n
        if check and n==0:
            a=s.copy();b=buf.copy();ref=np.zeros((20,26))
            nr,kr=P._tick(a,b,cmd,c,g,sc['q_field'],sc['E_field'],sc['ribs'],sc['h_s'],sc['xcmd'],sc['u'],sc['fext'],sc['act'],sc['rip'],alpha,tau,noise,kk,ref,0)
        n,kk=tick(s,buf,cmd,c,g,sc['q_field'],sc['E_field'],sc['ribs'],sc['h_s'],sc['xcmd'],sc['u'],sc['fext'],sc['act'],sc['rip'],alpha,tau,noise,kk,log,n)
        assert n>=0,'Plant failed'
        if check and n0==0:
            assert nr==n and kr==kk and np.array_equal(a,s) and np.array_equal(b,buf)
            assert np.array_equal(ref[:nr],log[:n,:26]),'Diagnostic tick changed physics'
        log[n0,5]=seen
    return sc,pl,start_state,c,command,log[:n],amplitude


def sense(sc,pl,ticks):
    pol=L.B.build('Huber-D4',L.O.chosen()['Huber-D4'],sc)
    pol.reset(dict(P.PLATFORM,tau_a=sc.get('tau_a_hat',pl['tau_a'])))
    sensing=pol.step if isinstance(pol.step,L.M.ForceSensing) else None
    calibrated=ticks.copy()
    def read(o):
        calibrated[o.n-1]=o.log[o.n-1]
        return 0.
    if sensing:sensing.step=read
    for n in range(1,len(ticks)+1):
        o=P.Obs(ticks,n)
        if sensing:sensing(o)
        else:read(o)
    return calibrated,pol


def estimator(sc,pl,cal,phi,without_h=False,learn=False):
    pol=L.B.build('Huber-D4',L.O.chosen()['Huber-D4'],sc)
    pol.reset(dict(P.PLATFORM,tau_a=sc.get('tau_a_hat',pl['tau_a'])))
    # Sensing has already run, so call the estimator directly.
    if without_h:pol.estimator.kind='L'
    pol.set_memory(phi)
    values=[];windows=[]
    for n in range(1,len(cal)+1):
        obs=SimpleNamespace(log=cal,n=n)
        kh,sig,chi=pol.estimator.estimate(obs,n,pol)
        if learn and np.isfinite(kh):pol.update_phase(obs,n)
        if cal[n-1,0]>=START and np.isfinite(kh):values.append(float(kh));windows.append([float(cal[n-1,0]),float(kh)])
    ll=cal if without_h else pol.estimator.private[:len(cal)]
    # The raw regression's weighted within-tick vectors, across the steady run.
    j0=max(pol.ndel+1,int(round(START/TC)));m=len(cal)-1-j0
    _,Y,Z,keep=L.common._build(ll,j0,m,0,pol.ns,pol.eb,pol.afc,pol.ndel,ll[-1,0],ll[-1,2],.192,True)
    t=cal[j0:j0+m,0]
    return dict(k_hat=summary(values),windows=windows,H_phi=float(pol.phi),
                H_tau=float(pol.tau_h),H_magnitude=float(pol.mag_h),
                weighted_Z_over_Y=describe(phasor(t[keep],Z[keep],CURRENT_FREQUENCY)/phasor(t[keep],Y[keep],CURRENT_FREQUENCY))),float(pol.phi)


CURRENT_FREQUENCY=20.8


def layers(data,frequency,phi):
    global CURRENT_FREQUENCY
    CURRENT_FREQUENCY=frequency
    sc,pl,state,c,command,hi,amplitude=data
    ticks=hi[::16,:26].copy();cal,_=sense(sc,pl,ticks)
    mask=hi[:,0]>=START;tm=hi[mask,0]
    def h(col):return phasor(tm,hi[mask,col],frequency)
    def ratio(num,den):return describe(h(num)/h(den))
    kc,ka,kt,op=constitutive(sc,state,c,frequency)
    _,_,k20,_=constitutive(sc,state,c,20.8)
    sampled=[]
    for offset in range(16):
        a=hi[offset::16];m=a[:,0]>=START;t=a[m,0]
        sampled.append(dict(offset_ms=offset*.5,**describe(phasor(t,a[m,18],frequency)/phasor(t,a[m,2],frequency))))
    q=ticks[:,0]>=START;t=ticks[q,0]
    cm=cal[:,0]>=START;tc=cal[cm,0]
    en_H,_=estimator(sc,pl,cal,phi)
    en_no,_=estimator(sc,pl,cal,0.,without_h=True)
    normal=h(3)/h(6);axial=h(18)/h(2)
    # The requested scalar closure uses the no-series, axial-projected K.
    ww=2*np.pi*frequency;mass=sc.get('tool_mass',0.);kr=sc['kr']
    requested=ka*(1+mass*ww*ww/kr)/(1+ka/kr)
    # Isolate extra static geometry and dynamic support at the level of this
    # small-signal formula; it is not an alternative simulation.
    theta=state[17];Zbody=sc['k']-sc['m']*ww*ww+1j*sc['c']*ww
    generalized=ka*(1+mass*ww*ww/kr)/(1+kc*op['series']+ka/Zbody)
    # Check the actual depth equation, including lateral arm compliance.
    geom=P._geometry(sc)
    surface=np.interp(hi[mask,26],sc['field_x0']+np.arange(len(sc['h_s']))*sc['dx'],sc['h_s'])
    lateral=hi[mask,19]*np.cos(hi[mask,17])-hi[mask,3]*np.sin(hi[mask,17])
    reconstructed=(hi[mask,2]+hi[mask,4]+surface-(hi[mask,18]+mass*hi[mask,27])/kr-c[26]*lateral)*np.cos(hi[mask,17])
    measured7=phasor(t,ticks[q,5],frequency)/phasor(t,command[q],frequency)
    true7=phasor(t,ticks[q,18],frequency)/phasor(t,command[q],frequency)
    comp7=phasor(tc,cal[cm,5],frequency)/phasor(tc,command[cm],frequency)
    result=dict(frequency=frequency,amplitude_command_m=amplitude,operating_point=op,
        L1_logged_col21=describe(k20),L1_logged_col21_measured=summary(hi[mask,21]),
        L1_formula_at_drive=describe(kt),L1_no_series_at_drive=describe(ka),
        L2_normal_force_over_depth=describe(normal),L2_constitutive_formula=describe(kc),
        L3_continuous_axial_over_encoder=describe(axial),L4_tick_axial_over_encoder=sampled[0],
        L4_sampler_offset=sampled,
        L5_compensated_force_over_aligned_encoder=describe(phasor(tc,cal[cm,5],frequency)/phasor(tc,cal[cm,2],frequency)),
        L5_compensated_force_over_physical_encoder=describe(phasor(tc,cal[cm,5],frequency)/phasor(t,ticks[q,2],frequency)),
        L6_estimator_H=en_H,L6_estimator_no_H=en_no,
        L7_command_to_raw_sensor=describe(measured7),L7_command_to_compensated_force=describe(comp7),
        L7_command_to_true_axial=describe(true7),encoder_over_command=describe(phasor(t,ticks[q,2],frequency)/phasor(t,command[q],frequency)),
        scalar_quasistatic_prediction=describe(requested),scalar_relative_closure_error=describe(axial/requested),
        geometry_support_prediction=describe(generalized),geometry_support_relative_closure_error=describe(axial/generalized),
        validation=dict(depth_equation_max_abs_error_m=float(np.max(abs(reconstructed-hi[mask,6]))),
            lateral_contact_excursion_m=float(np.ptp(hi[mask,26])),
            relative_K_span=float(np.ptp(hi[mask,21])/np.mean(hi[mask,21])),
            force_dc=float(np.mean(hi[mask,18])),force_amplitude=float(abs(h(18))),
            memory_initial_max=float(np.max(abs(state[37:45])))))
    return result


def image_analysis(data,frequency=20.8):
    sc,pl,state,c,command,hi,amplitude=data
    q=hi[:,0]>=START;t=hi[q,0];m=sc.get('tool_mass',0.);kr=sc['kr']
    ysum=0j;fqsum=0j;fmsum=0j;factual=0j;contributions=[]
    fdqsum=0j;fdmsum=0j
    _,_,_,op=constitutive(sc,state,c,frequency)
    alpha,tau,_=P.memory_arrays(sc);decay=np.exp(-DT/tau)
    # All 16 distinct aliases below the 2 kHz logging Nyquist. Fit jointly to
    # avoid cross-leakage of the real signal's conjugates at fractional cycles.
    indices=list(range(-8,8));frequencies=[frequency+n/TC for n in indices]
    s=t-t.mean();design=np.column_stack([v for nu in frequencies for v in (np.cos(2*np.pi*nu*t),np.sin(2*np.pi*nu*t))]+[np.ones(len(t)),s,s*s])
    coef=np.linalg.lstsq(design,np.c_[hi[q,2],hi[q,18]],rcond=None)[0]
    for j,n in enumerate(indices):
        nu=frequency+n/TC
        Y=complex(coef[2*j,0],-coef[2*j+1,0]);Factual=complex(coef[2*j,1],-coef[2*j+1,1])
        _,K,_,_=constitutive(sc,state,c,abs(nu))
        if nu<0:K=K.conjugate()
        ww=2*np.pi*nu
        Xq=(1+m*ww*ww/kr)/(1+K/kr)
        Xm=kr/(kr+K-m*ww*ww)
        ysum+=Y;fqsum+=K*Xq*Y;fmsum+=K*Xm*Y;factual+=Factual
        # Same alternatives with the actual 0.5 ms difference operators and
        # exact implemented Prony ramp weights. This avoids extrapolating
        # continuous omega^2 up to the logging Nyquist in the alias sum.
        z=np.exp(-1j*ww*DT);omega2=4*np.sin(ww*DT/2)**2/DT**2
        Kd=op['normal_elastic_tangent']*(1+np.sum(alpha*tau/DT*(1-decay)*(1-z)/(1-decay*z)))*np.cos(state[17])**2
        body=sc['k']-sc['m']*omega2+sc['c']*(1-z)/DT
        effective=Kd/(1+Kd/body)
        Gdq=effective*(1+m*omega2/kr)/(1+effective/kr)
        Gdm=effective*kr/(kr+effective-m*omega2)
        fdqsum+=Gdq*Y;fdmsum+=Gdm*Y
        contributions.append(dict(n=n,frequency=nu,Y=describe(Y),actual_F=describe(Factual),
                                  quasistatic_x_over_y=describe(Xq),mass_spring_x_over_y=describe(Xm),
                                  quasi_F=describe(K*Xq*Y),spring_F=describe(K*Xm*Y),
                                  discrete_quasi_F=describe(Gdq*Y),discrete_spring_F=describe(Gdm*Y)))
    rows=[]
    for f in (104.,125.-frequency):
        _,K,_,_=constitutive(sc,state,c,f);w=2*np.pi*f
        rows.append(dict(frequency=f,K_no_series=describe(K),
            quasi_x_over_y=describe((1+m*w*w/kr)/(1+K/kr)),
            mass_spring_x_over_y=describe(kr/(kr+K-m*w*w))))
    return dict(natural_frequency_Hz=math.sqrt(kr/m)/(2*np.pi) if m else None,
        at_104=rows,harmonics=contributions,
        predicted_sampled_quasistatic_scalar=describe(fqsum/ysum),
        predicted_sampled_true_mass_spring_scalar=describe(fmsum/ysum),
        discrete_1d_quasistatic=describe(fdqsum/ysum),
        discrete_1d_mass_spring=describe(fdmsum/ysum),
        reconstructed_measured_sampled=describe(factual/ysum),
        note='Derived 1D arm alternatives at the same measured encoder harmonics, not a changed-plant simulation. Continuous scalar sums omit geometry/support/discretization; discrete_1d adds implemented Prony, acceleration differences, and body support, but omits lateral coupling/spatial gradients. Its measured quasi closure bounds that residual before interpreting the mass-spring alternative.')


def finish(result):
    for condition in result['conditions']:
        for channel in ('L7_command_to_raw_sensor','L7_command_to_compensated_force','L7_command_to_true_axial','encoder_over_command'):
            angles=np.unwrap(np.array([r[channel]['phase_deg'] for r in condition['sweep']])*np.pi/180)*180/np.pi
            for r,angle in zip(condition['sweep'],angles):r[channel]['unwrapped_phase_deg']=float(angle)
            base=next(r for r in condition['sweep'] if r['frequency']==20.8)[channel]['unwrapped_phase_deg']
            for r in condition['interventions']:
                angle=r[channel]['phase_deg'];r[channel]['unwrapped_phase_deg']=float(angle+360*round((base-angle)/360))
        for r in condition['sweep']+condition['interventions']:
            k0=r['L1_logged_col21']['gain']
            r['gain_ratios_to_logged_K20']={k:r[k]['gain']/k0 for k in (
                'L1_no_series_at_drive','L2_normal_force_over_depth','L3_continuous_axial_over_encoder',
                'L4_tick_axial_over_encoder','L5_compensated_force_over_aligned_encoder','L7_command_to_compensated_force')}
            r['gain_ratios_to_logged_K20'].update(k_hat_H=r['L6_estimator_H']['k_hat']['median']/k0,
                                                  k_hat_no_H=r['L6_estimator_no_H']['k_hat']['median']/k0)
        base=next(r for r in condition['sweep'] if r['frequency']==20.8)
        ref=base['L1_logged_col21']['gain']
        condition['images']['sampled_scalar_quasi_ratio_to_K20']=condition['images']['predicted_sampled_quasistatic_scalar']['gain']/ref
        condition['images']['sampled_scalar_mass_spring_ratio_to_K20']=condition['images']['predicted_sampled_true_mass_spring_scalar']['gain']/ref
        if 'discrete_1d_quasistatic' in condition['images']:
            condition['images']['discrete_quasi_ratio_to_K20']=condition['images']['discrete_1d_quasistatic']['gain']/ref
            condition['images']['discrete_mass_spring_ratio_to_K20']=condition['images']['discrete_1d_mass_spring']['gain']/ref
    result['metadata']['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def main():
    assert not OUT.exists(),'Refuse to overwrite an existing result'
    begin=time.perf_counter();source=json.loads(SOURCE.read_text(encoding='utf-8'));rows=[]
    for label in LABELS:
        row=next(r for r in source['episodes'] if r['label']==label)
        sc,origin=scene(row)
        base20=run(sc,20.8,check=True)
        cal,_=sense(base20[0],base20[1],base20[5][::16,:26])
        # H's original synchronous phase calculation at its design frequency;
        # no noise-precision gate is imposed during this noiseless calibration.
        global CURRENT_FREQUENCY
        CURRENT_FREQUENCY=20.8
        calibration,phi=estimator(sc,base20[1],cal,0.,without_h=True,learn=True)
        swept=[]
        for f in FREQUENCIES:
            data=base20 if f==20.8 else run(sc,f)
            swept.append(layers(data,f,phi))
        interventions=[]
        amplitude=base20[-1]
        for change in ('mass0','kr_large','fast1ms','fast4ms','offset0'):
            data=run(sc,20.8,change,amplitude=amplitude)
            interventions.append(dict(change=change,**layers(data,20.8,phi)))
        phases=[]
        for phase in (0.,np.pi/2,np.pi,3*np.pi/2):
            data=base20 if phase==0 else run(sc,20.8,phase=phase,amplitude=amplitude)
            log=data[5][::16];q=log[:,0]>=START;t=log[q,0]
            phases.append(dict(input_phase_deg=float(phase*180/np.pi),
                **describe(phasor(t,log[q,18],20.8)/phasor(t,log[q,2],20.8))))
        # Half-amplitude check bounds local nonlinearity in the main diagnostic.
        half=run(sc,20.8,amplitude=amplitude/2)
        hh=half[5][::16];q=hh[:,0]>=START
        linearity=describe(phasor(hh[q,0],hh[q,18],20.8)/phasor(hh[q,0],hh[q,2],20.8))
        rows.append(dict(label=label,source_case=row['case'],origin=origin,
            parameters={k:L.clean(sc.get(k)) for k in ('tool_mass','tool_mass_hat','kr','rho','K_t','m','k','c','alpha','tau','sensor_delay_ms','sensor_delay_hat','tau_a_true','tau_a_hat','phi_hc')},
            H_calibration_phase=phi,H_calibration=calibration,sweep=swept,
            interventions=interventions,input_phase_sweep=phases,
            half_amplitude_L4=linearity,images=image_analysis(base20)))
        print('completed',label,round(time.perf_counter()-begin,1),flush=True)
    result=dict(metadata=dict(noise=False,open_loop=True,highrate_dt=DT,control_dt=TC,
        duration=DURATION,fit_start=START,frequencies=FREQUENCIES,
        operating_point='Original t=4s scan command held constant. Original spatial material map unchanged and static; no scan, breathing, activation, external load or ripple. Memory starts exactly relaxed. Small lateral excursions and K modulation measured.',
        H='Calibrated once per condition by original update_phase on noise-free no-H20.8Hz data, then held fixed for all frequencies/interventions; no precision gate in this diagnostic calibration. H increment default; no gates or memory filters.',
        phase='output relative to input; scalar k_hat has no phase. weighted_Z_over_Y phase is supplied separately, not assigned to scalar k_hat.',
        interventions='Physical fast pole only, keeping total tau_a fixed; ForceSensing retains nominal2ms and original calibration. mass0 removes true and compensated tool mass; offset0 removes true and calibrated offset. kr_large=1e12.',
        L1='col21/22 always20.8Hz. Reevaluate original constitutive/series formula at every drive frequency as a separate field.',
        L7='Three explicit channels: command to raw sensed force, ForceSensing-compensated force consumed by Method, and true axial force. Command is the held sequence supplied to plant.',
        plant_source_sha256=hashlib.sha256(Path(P.__file__).read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),conditions=rows,seconds=time.perf_counter()-begin)
    finish(result)
    OUT.write_text(json.dumps(L.clean(result),ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')
    print('saved',OUT,'seconds',result['seconds'],flush=True)


if __name__=='__main__':main()
