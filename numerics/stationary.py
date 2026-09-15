"""Matched-precision stationary responses and fixed-time gain partials."""
from __future__ import annotations
import csv,hashlib,json
from pathlib import Path
import mpmath as mp
from .fixed_model import fixed_model,INPUTS
GAINS=(10, 25, 30, 50, 75, 100, 125, 140, 150, 160, 175, 200, 300, 400, 600, 1000, 3000, 10000, 100000, 1000000, 10000000, 100000000)
METRICS=('normalized','raw')
DIGITS=(80,110)
CORRECTION_FLOOR='1e-8'
def require(test: bool, message: str) -> None:
    if not test: raise ValueError(message)

def top_pair(T):
    values,V=mp.eigsy(T.T*T)
    require(values[2]>values[1] and values[2]>0,'Top singular value not simple.')
    G=mp.sqrt(values[2]);v=V[:,2];u=T*v/G
    return G,u,v,values,V

def discrepancy_values(G,t,D,Ghat,that,Dhat,gamma,r,L,d):
    """Percent errors and signed residuals with explicit omission rules.
    Omit correction-relative error when |gamma-Ghat| <= 1e-8 gamma.
    Omit derivative-relative error at D=0 and scaled residuals at L=0.
    """
    require(G>0 and t>0 and gamma>0 and r>0 and d>0,'Invalid response inputs.')
    eG=G-Ghat;et=t-that;eD=D-Dhat;corr=abs(gamma-Ghat)
    corr_ok=corr>mp.mpf(CORRECTION_FLOOR)*gamma
    return dict(height_error_pct=100*abs(eG)/G,
      correction_error_pct=100*abs(eG)/corr if corr_ok else None,
      time_error_pct=100*abs(et)/abs(t),
      sensitivity_error_pct=100*abs(eD)/abs(D) if D!=0 else None,
      projection_error_pct=100*abs(G-gamma)/G,
      height_residual=eG,time_residual=et,sensitivity_residual=eD,
      R_G=eG*r*r/(L*L) if L!=0 else None,
      R_t=et*r*r/L if L!=0 else None,
      R_D=eD*r**3/(d*L*L) if L!=0 else None,
      correction_denominator=corr,correction_ratio_to_gamma=corr/gamma,
      correction_defined=corr_ok,derivative_relative_defined=bool(D!=0))

def evaluate(metric: str, gain: int, digits: int):
    dat=fixed_model(metric);k=mp.mpf(gain);r=dat['d']*k
    M=dat['A0']-r*dat['Q'];N=M/r;L=mp.log(dat['chi']*r/dat['a'])
    require(L>0,'Selected root seeds require positive L.')
    cache={}
    def state(tau):
        key=mp.nstr(tau,mp.mp.dps)
        if key not in cache:
            T=mp.expm(N*tau);G,u,v,values,V=top_pair(T)
            cache[key]=(T,G,u,v,values,V,(u.T*N*T*v)[0])
        return cache[key]
    lo=max(mp.mpf('0.0001'),L-mp.mpf('0.5'));hi=L+mp.mpf('0.5')
    for _ in range(18):
        if state(lo)[-1]>0 and state(hi)[-1]<0:break
        lo/=2;hi+=1
    else:raise ValueError('No descending stationary bracket located.')
    tau=mp.findroot(lambda z:state(z)[-1],(lo,hi),solver='anderson',
        tol=mp.mpf(10)**(-digits+15),maxsteps=100,verify=True)
    require(lo<tau<hi and tau>0,'Root left positive bracket.')
    T,G,u,v,values,V,stationarity=state(tau);t=tau/r
    # Exact simple-eigenvalue differentiation of T^T T for time curvature.
    Tt=M*T;Ttt=M*Tt;Xt=Tt.T*T+T.T*Tt
    Xtt=Ttt.T*T+2*Tt.T*Tt+T.T*Ttt
    lam1=(v.T*Xt*v)[0]
    lam2=(v.T*Xtt*v)[0]+2*sum((V[:,j].T*Xt*v)[0]**2/(values[2]-values[j]) for j in range(2))
    curvature=lam2/(2*G)-lam1**2/(4*G**3)
    block=mp.zeros(6)
    for i in range(3):
        for j in range(3):
            block[i,j]=t*M[i,j];block[i+3,j+3]=t*M[i,j]
            block[i,j+3]=-t*dat['d']*dat['Q'][i,j]
    Tk=mp.expm(block)[:3,3:];D=(u.T*Tk*v)[0]
    Ghat=dat['gamma']-(dat['a']*(L+1)-dat['h'])/r
    that=L/r;Dhat=dat['d']*(dat['a']*L-dat['h'])/r**2
    errs=discrepancy_values(G,t,D,Ghat,that,Dhat,dat['gamma'],r,L,dat['d'])
    residual=max(abs(stationarity),mp.norm(T*v-G*u),mp.norm(T.T*u-G*v))
    require(residual<mp.mpf('1e-60'),'Stationary or singular-pair residual failed.')
    require(curvature<0,'Nonnegative stationary curvature.')
    gap=mp.sqrt(values[2])-mp.sqrt(max(mp.mpf(0),values[1]))
    step=mp.mpf('0.001')
    require(state(tau-step)[-1]>0 and state(tau+step)[-1]<0,'Local sign check failed.')
    # Independent 3x3 spectral exponential and divided-difference Frechet
    # calculation check the primary 3x3/6x6 direct matrix exponentials.
    eigenvalues,W=mp.eig(M);Wi=W**-1
    Teig=W*mp.diag([mp.exp(q*t) for q in eigenvalues])*Wi
    Z=Wi*(-dat['d']*dat['Q'])*W;F=mp.zeros(3)
    for i in range(3):
        for j in range(3):
            if i==j:div=t*mp.exp(eigenvalues[i]*t)
            else:
                den=eigenvalues[i]-eigenvalues[j]
                require(abs(den)>mp.mpf('1e-20'),'Coalescing spectral eigenvalues.')
                div=mp.exp(eigenvalues[j]*t)*mp.expm1(den*t)/den
            F[i,j]=Z[i,j]*div
    spectral_D=(u.T*W*F*Wi*v)[0]
    propagation_error=mp.norm(Teig-T)/max(1,mp.norm(T))
    sensitivity_error=abs(spectral_D-D)/abs(D)
    require(propagation_error<mp.mpf('1e-55') and sensitivity_error<mp.mpf('1e-55'),
            'Spectral/block agreement failed.')
    show=lambda z:mp.nstr(z,digits-10)
    numeric=dict(r=r,L=L,peak=G,peak_time=t,gain_derivative=D,predicted_peak=Ghat,
       predicted_time=that,predicted_gain_derivative=Dhat,gamma=dat['gamma'],a_p=dat['a'],
       h_p=dat['h'],chi=dat['chi'],j_p=dat['j'],d=dat['d'],stationarity_residual=residual,
       curvature=curvature,singular_gap=gap,bracket_lower=lo,bracket_upper=hi,
       spectral_propagator_relative_error=propagation_error,
       spectral_sensitivity_relative_error=sensitivity_error)
    row=dict(metric=metric,gain=gain,dps=digits,
        status='stationary witness; not an interval-certified global extremum',
        log_parameter_below_one=bool(L<1),
        below_supplied_comparison_range=(bool(gain<INPUTS['comparison_inputs']['rate']) if metric=='normalized' else None),
        sensitivity_opposite_sign=bool(D*Dhat<0))
    row.update({key:show(val) for key,val in numeric.items()})
    row.update({key:(val if isinstance(val,bool) or val is None else show(val)) for key,val in errs.items()})
    return row

def run(output: Path):
    runs=[]
    for digits in DIGITS:
        mp.mp.dps=digits;rows=[]
        for metric in METRICS:
            for gain in GAINS:
                rows.append(evaluate(metric,gain,digits))
                print(f'{digits} digits {metric} k={gain}: stationary and spectral checks pass',flush=True)
        runs.append(dict(dps=digits,rows=rows))
    mp.mp.dps=DIGITS[-1];comparisons=[]
    for lo,hi in zip(runs[0]['rows'],runs[1]['rows']):
        for key in ('peak','peak_time','gain_derivative','R_G','R_t','R_D'):
            den=abs(mp.mpf(hi[key]));require(den>0,'Zero precision denominator.')
            err=abs(mp.mpf(lo[key])-mp.mpf(hi[key]))/den
            require(err<mp.mpf('1e-55'),f'Matched precision failed for {key}.')
            comparisons.append(dict(metric=hi['metric'],gain=hi['gain'],quantity=key,
                                    relative_difference=mp.nstr(err,25)))
    constants={}
    for metric in METRICS:
        D=fixed_model(metric);a=D['a'];h=D['h'];d=D['d'];chi=D['chi']
        constants[metric]={key:mp.nstr(D[key],90) for key in ('gamma','chi','a','h','j','d')}
        constants[metric].update({
            'gain_L_equals_1':mp.nstr(mp.e*a/(d*chi),90),
            'gain_L_equals_2':mp.nstr(mp.e**2*a/(d*chi),90),
            'gain_leading_sensitivity_zero':mp.nstr(a*mp.exp(h/a)/(d*chi),90),
            'gain_predicted_height_correction_zero':mp.nstr(a*mp.exp(h/a-1)/(d*chi),90)})
    obj=dict(status='passed',scope='numerical stationary diagnostics only; no new interval global certificate',
      methods=dict(precisions=list(DIGITS),sensitivity='fixed-time 6x6 block Frechet exponential',
         independent_crosscheck='3x3 spectral exponential and exponential divided-difference Frechet formula',
         root='positive descending bracket in fast time; negative curvature and local sign checks',
         correction_floor_relative_to_gamma=CORRECTION_FLOOR,matched_relative_tolerance='1e-55'),
      gains=list(GAINS),metrics=list(METRICS),constants=constants,runs=runs,
      precision_comparisons=comparisons,stationary_evaluations=sum(len(x['rows']) for x in runs),
      matched_precision_comparisons=len(comparisons),
      independent_spectral_checks=2*sum(len(x['rows']) for x in runs),
      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(obj,indent=2)+'\n')
    with output.with_suffix('.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(runs[-1]['rows'][0]));w.writeheader();w.writerows(runs[-1]['rows'])
    print('PASS:',obj['stationary_evaluations'],'stationary evaluations;',len(comparisons),'precision comparisons.',flush=True)
    return obj
