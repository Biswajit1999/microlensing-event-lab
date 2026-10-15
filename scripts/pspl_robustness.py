"""Additional diagnostics for the OGLE PSPL fit; no lens-mass inference.

Requires the repository's existing fit_ogle_event.py and numpy/scipy. This
module never changes the browser's baseline fit JSON unless explicitly run.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from fit_ogle_event import REFERENCE_MAG, SOURCE, fit_model, flux_model


PARAMETER_NAMES = ('t0_hjd','tE_days','u0','source_flux','blend_flux')


def analytic_flux_jacobian(time, parameters):
    """Exact derivative dF/d(t0,tE,u0,Fs,Fb) for point-source/point-lens.

    The Jacobian is in *physical* parameter units. u0 and tE must be >0.
    """
    t0, te, u0, fs, fb = map(float, parameters)
    if te <= 0 or u0 <= 0 or fs < 0 or fb < 0:
        raise ValueError('invalid physical PSPL parameters')
    t = np.asarray(time, dtype=float)
    tau = (t-t0)/te
    u = np.hypot(u0,tau)
    a = (u*u+2)/(u*np.sqrt(u*u+4))
    da_du = -8/(u*u*(u*u+4)**1.5)
    d_t0 = fs*da_du*(-tau/(u*te))
    d_te = fs*da_du*(-tau*tau/(u*te))
    d_u0 = fs*da_du*(u0/u)
    return np.column_stack((d_t0,d_te,d_u0,a,np.ones_like(t)))


def fisher_covariance(time, error, parameters, *, rcond=1e-10):
    """Local Fisher covariance, assuming independent Gaussian flux errors.

    SVD on column-scaled design matrix guards against spurious rank estimates
    caused by different physical parameter units. Rank-deficient fits do not
    receive a misleading finite uncertainty estimate.
    """
    sigma = np.asarray(error, dtype=float)
    if np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError('positive finite errors required')
    j = analytic_flux_jacobian(time, parameters)/sigma[:,None]
    scale = np.linalg.norm(j,axis=0)
    if np.any(scale <= 0):
        return {'identifiable':False,'rank':int(np.count_nonzero(scale)), 'formal_sigma':None}
    _,singular,vt = np.linalg.svd(j/scale,full_matrices=False)
    rank = int(np.sum(singular > singular[0]*rcond))
    if rank < len(PARAMETER_NAMES):
        return {'identifiable':False,'rank':rank,'formal_sigma':None}
    inverse_scaled = (vt.T/(singular**2))@vt
    covariance = inverse_scaled/scale[:,None]/scale[None,:]
    return {'identifiable':True,'rank':rank,
            'scaled_condition_number':float(singular[0]/singular[-1]),
            'formal_sigma':dict(zip(PARAMETER_NAMES,np.sqrt(np.diag(covariance)).tolist())),
            'covariance':covariance.tolist()}


def observing_condition_diagnostics(flux, model, error, seeing, sky):
    """Noncausal seeing/sky association and high-seeing residual scatter.

    Correlations are exploratory and ignore temporal covariance/multiple tests.
    """
    f, m, s, seeing, sky = map(lambda x:np.asarray(x,dtype=float),
                                (flux,model,error,seeing,sky))
    if not (len(f)==len(m)==len(s)==len(seeing)==len(sky)) or len(f)<8:
        raise ValueError('at least eight matched observations required')
    if np.any(s<=0) or not all(np.all(np.isfinite(x)) for x in (f,m,s,seeing,sky)):
        raise ValueError('invalid observations')
    z=(f-m)/s
    q=float(np.quantile(seeing,0.75))
    high=seeing>q
    if not (np.any(high) and np.any(~high)):
        raise ValueError('insufficient seeing variation')
    def corr(x):
        v=spearmanr(np.abs(z),x).statistic
        return float(v) if np.isfinite(v) else None
    return {'n':len(z),'high_seeing_threshold_px':q,
            'high_seeing_count':int(high.sum()),
            'high_seeing_mean_squared_z':float(np.mean(z[high]**2)),
            'other_mean_squared_z':float(np.mean(z[~high]**2)),
            'spearman_abs_z_seeing':corr(seeing),
            'spearman_abs_z_sky':corr(sky)}


def parametric_bootstrap(time, error, parameters, *, n_draws=64, seed=1029):
    """Refit Gaussian synthetic flux around fitted PSPL curve.

    This estimates conditional noise-only scatter, not uncertainty due to
    correlated systematics, blending priors, or alternate lens models.
    """
    if n_draws < 8:
        raise ValueError('at least eight draws required')
    t=np.asarray(time,dtype=float)
    s=np.asarray(error,dtype=float)
    if np.any(s<=0) or np.any(~np.isfinite(s)):
        raise ValueError('positive finite errors required')
    baseline=flux_model(t,*parameters)
    rng=np.random.default_rng(seed)
    samples=[]
    failures=0
    for _ in range(n_draws):
        noisy=baseline+rng.normal(0,s)
        try:
            fit=fit_model(t,noisy,s,blended=True)
            p=fit['parameters']
            samples.append([p[name] for name in PARAMETER_NAMES])
        except (ValueError,RuntimeError,FloatingPointError):
            failures+=1
    if len(samples)<8:
        raise RuntimeError('too few successful refits')
    values=np.asarray(samples)
    return {'seed':seed,'draws_requested':n_draws,'draws_successful':len(samples),
            'failures':failures,
            'central_68_percentile_intervals':{
                name:[float(v) for v in np.quantile(values[:,i],[0.16,0.84])]
                for i,name in enumerate(PARAMETER_NAMES)}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=SOURCE)
    parser.add_argument('--output',type=Path,default=SOURCE.parent/'ogle-2023-blg-0100-robustness.json')
    parser.add_argument('--bootstrap',type=int,default=64)
    args=parser.parse_args()
    raw=np.loadtxt(args.input)
    t,mag,magerr,seeing,sky=raw.T
    f=10**(-0.4*(mag-REFERENCE_MAG))
    s=f*np.log(10)*0.4*magerr
    fit=fit_model(t,f,s,blended=True)
    p=tuple(fit['parameters'][key] for key in PARAMETER_NAMES)
    output={'event':'OGLE-2023-BLG-0100','n_observations':len(t),
            'limitations':'Independent Gaussian errors; PSPL model conditional; no unique lens mass.',
            'fisher':fisher_covariance(t,s,p),
            'conditions':observing_condition_diagnostics(f,fit['prediction'],s,seeing,sky),
            'bootstrap':parametric_bootstrap(t,s,p,n_draws=args.bootstrap)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n',encoding='utf8')
    print(f'Wrote PSPL diagnostics for {len(t)} observations to {args.output}')


if __name__=='__main__':
    main()
