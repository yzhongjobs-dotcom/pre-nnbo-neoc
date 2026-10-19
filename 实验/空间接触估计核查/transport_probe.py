"""Conditional force/moment transport inversion, with explicit model mismatch.

This is an identifiability/mechanism probe, not a task controller or physical
validation. It refines B0's existing joint material/surface estimation idea.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent
A,L,R=.01,.02,.05


def moments(d):
    a=min(A,np.sqrt(max(0.,2*R*d)))
    aa=L*(2*a*d-a**3/(3*R))
    bb=L*(2*a**3*d/3-a**5/(5*R))
    return aa,bb


def inverse_load(f,k):
    return brentq(lambda d:k*moments(d)[0]-f,1e-7,.1,xtol=1e-14)


def main():
    dt=.01;times=np.arange(0, .064/.012,dt);speed=.012
    nodes,weights=leggauss(257);u=A*nodes;weights=L*A*weights
    field_cases=[('affine_softening',-3e7,False),('affine_hardening',3e7,False),
                 ('sharp_transition',0.,False),('affine_with_viscosity',3e7,True)]
    rows=[];arrays={'times':times}
    for index,(name,gradient,viscous) in enumerate(field_cases):
        def field(x):
            return 4.5e6*(1+.4*np.tanh((x-.032)/.004)) if name=='sharp_transition' else 4.5e6+gradient*x
        height=.005*np.sin(2*np.pi*.2*times)
        true_k=field(speed*times)
        indentation=[]
        for t in times:
            k=field(speed*t+u)
            indentation.append(brentq(lambda d:np.dot(weights,k*np.maximum(0.,d-u*u/(2*R)))-3.,1e-7,.03,xtol=1e-13))
        indentation=np.array(indentation);robot=height+indentation
        rate=np.gradient(indentation,dt);forces=[];torques=[]
        for j,t in enumerate(times):
            delta=indentation[j]-u*u/(2*R);k=field(speed*t+u)
            p=k*delta
            if viscous:p+=k*.008*(rate[j]+u*speed/R)
            p=np.where(delta>0,np.maximum(p,0.),0.)
            forces.append(np.dot(weights,p));torques.append(np.dot(weights,u*p))
        forces=np.array(forces);torques=np.array(torques)
        for noise_mode in ['clean','corrected_moment_sensitivity','raw_wrench']:
            rng=np.random.default_rng(9401+index)
            f=forces.copy();mm=torques.copy()
            if noise_mode!='clean':f+=rng.normal(0,.03,len(f))
            if noise_mode=='corrected_moment_sensitivity':mm+=rng.normal(0,.0001,len(f))
            if noise_mode=='raw_wrench':
                # Same raw force/torque noise convention as the control carrier.
                ft=.2*forces;raw_my=-torques+.025*ft
                measured_ft=ft+rng.normal(0,.03,len(f))
                measured_my=raw_my+rng.normal(0,.0001,len(f))
                mm=-measured_my+.025*measured_ft
            kh=np.empty(len(times));hh=np.empty(len(times));fixed=np.empty(len(times))
            # Same declared initial surface reference as B0; no hidden later height.
            kh[0]=f[0]/moments(robot[0])[0]
            clips=0
            for j in range(len(times)):
                d=inverse_load(max(.05,f[j]),kh[j]);hh[j]=robot[j]-d
                fixed[j]=robot[j]-inverse_load(max(.05,f[j]),kh[0])
                if j+1<len(times):
                    derivative=speed*mm[j]/max(moments(d)[1],1e-16)
                    updated=kh[j]+dt*derivative
                    kh[j+1]=np.clip(updated,1e5,2e7);clips+=int(updated!=kh[j+1])
            key=name+'_'+noise_mode
            arrays[key]=np.column_stack([true_k,kh,height,hh,fixed,forces,torques,robot])
            row={'case':name,'noise_mode':noise_mode,'relative_material_RMS':float(np.sqrt(np.mean(((kh-true_k)/true_k)**2))),
                 'height_RMS_mm':float(np.sqrt(np.mean((hh-height)**2))*1000),
                 'fixed_material_height_RMS_mm':float(np.sqrt(np.mean((fixed-height)**2))*1000),
                 'parameter_projection_samples':clips}
            rows.append(row);print(row,flush=True)
    out={'rows':rows,'sample_Hz':100,'motion_amplitude_mm':5,'motion_frequency_Hz':.2,
         'noise_std_Fn_Ft_N':.03,'noise_std_raw_My_Nm':.0001,
         'corrected_normal_moment_std_Nm':float(np.sqrt(.0001**2+(.025*.03)**2)),
         'initialization':'Initial surface height and force/robot position; no hidden later material or motion in inversion; declared known tangential moment arm',
         'timing':'Aligned sample-time estimates. A 20ms sensor pipeline would make them 20ms old; no current-state predictor/performance claim.',
         'scope':'Elastic affine field is the exact-model identification probe. Sharp boundary and convective viscosity deliberately violate it. Same physical-law data cannot validate real tissue; no task-level or novelty claim.',
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'transport_probe.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(HERE/'transport_probe.npz',**arrays)


if __name__=='__main__':main()
