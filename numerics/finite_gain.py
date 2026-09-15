"""Finite sampled gain searches, stationary candidates and direct spectra.
Results are numerical candidates rather than interval-certified optimizers.
"""
from __future__ import annotations
import csv,hashlib,json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm,eigvals,svdvals
from scipy.optimize import minimize_scalar,brentq
STATUS='numerical stationary witness; not an interval-certified global extremum'
def require(ok:bool,msg:str)->None:
    if not ok:raise ValueError(msg)

def matrix_float(k:float):
    require(np.isfinite(k) and k>0,'Gain must be positive and finite.')
    return np.array([[0.,5.,0.],[0.,0.,0.],[.5,0.,-5.]])+k*np.outer([1.,.6,0.],[-1.,0.,1.])

def grid_peak(k:float,density:int=1,horizon:float=24.):
    """Refine all sampled humps and check sampled semigroup tail reduction.
    This is a floating-point diagnostic; none of the time cells is enclosed.
    """
    require(density in (1,2),'Supported grid density is 1 or 2.')
    require(horizon>0,'Horizon must be positive.')
    A=matrix_float(k)
    ts=np.unique(np.r_[0,np.geomspace(1e-6,horizon,180*density),np.linspace(0,horizon,140*density)])
    Es=expm(ts[:,None,None]*A)
    gs=np.linalg.svd(Es,compute_uv=False)[:,0]
    require(float(gs[-1])<1,'Sampled tail norm not below one.')
    maxima=[j for j in range(1,len(ts)-1) if gs[j]>=gs[j-1] and gs[j]>=gs[j+1]]
    candidates=[(1.,0.)]
    def response(t):
        E=expm(t*A);u,w,vh=np.linalg.svd(E)
        return float(w[0]),float(u[:,0]@(A@E)@vh[0,:])
    for j in maxima:
        lo,hi=float(ts[j-1]),float(ts[j+1])
        if response(lo)[1]>0 and response(hi)[1]<0:
            t=brentq(lambda t:response(t)[1],lo,hi,xtol=5e-15)
        else:
            res=minimize_scalar(lambda t:-response(t)[0],bounds=(lo,hi),method='bounded',options={'xatol':1e-13})
            require(res.success,'Time refinement failed.');t=float(res.x)
        candidates.append((response(t)[0],t))
    G,t=max(candidates)
    require(t>0,'No positive-time peak found.')
    return dict(gain=float(k),peak=float(G),time=float(t),alpha=float(max(eigvals(A).real)),
        sample_count=int(len(ts)),sampled_humps=len(maxima),sampled_tail_norm=float(gs[-1]),
        horizon=float(horizon),density=density,status=STATUS)

def mp_matrix(k):
    B=mp.matrix([-1,-mp.mpf(3)/5,0])*mp.matrix([1,0,-1]).T
    M=mp.matrix([[0,5,0],[0,0,0],[mp.mpf(1)/2,0,-5]])+k*B
    return M,B

def mp_state(k,t,derivatives=True,second=False):
    M,B=mp_matrix(k);T=mp.expm(t*M);values,V=mp.eigsy(T.T*T)
    G=mp.sqrt(values[2]);v=V[:,2];u=T*v/G;Tt=M*T
    result=dict(G=G,u=u,v=v,T=T,Gt=(u.T*Tt*v)[0],gap=G-mp.sqrt(values[1]))
    if not derivatives:return result
    block=mp.zeros(9 if second else 6)
    for j in range(3 if second else 2):
        block[3*j:3*j+3,3*j:3*j+3]=t*M
        if j<(2 if second else 1):block[3*j:3*j+3,3*j+3:3*j+6]=t*B
    Eb=mp.expm(block);Tk=Eb[:3,3:6];Gk=(u.T*Tk*v)[0]
    Xt=Tt.T*T+T.T*Tt;Xtt=(M*Tt).T*T+2*Tt.T*Tt+T.T*(M*Tt)
    def curvature(Xi,Xj,Xij,li,lj):
        val=(v.T*Xij*v)[0]+2*sum((V[:,h].T*Xi*v)[0]*(V[:,h].T*Xj*v)[0]/(values[2]-values[h]) for h in range(2))
        return val/(2*G)-li*lj/(4*G**3)
    lt=(v.T*Xt*v)[0];Gtt=curvature(Xt,Xt,Xtt,lt,lt)
    result.update(Gk=Gk,Gtt=Gtt,Tk=Tk)
    if second:
        Tkk=2*Eb[:3,6:9];Ttk=B*T+M*Tk;Xk=Tk.T*T+T.T*Tk
        Xkk=Tkk.T*T+2*Tk.T*Tk+T.T*Tkk
        Xtk=Ttk.T*T+Tt.T*Tk+Tk.T*Tt+T.T*Ttk
        lk=(v.T*Xk*v)[0]
        Gkk=curvature(Xk,Xk,Xkk,lk,lk);Gtk=curvature(Xt,Xk,Xtk,lt,lk)
        result.update(Gkk=Gkk,Gtk=Gtk,branch_second_derivative=Gkk-Gtk**2/Gtt)
    return result

def high_precision_row(name,k,dps,tseed=None,second=False):
    if tseed is None:tseed=mp.mpf(str(grid_peak(float(k))['time']))
    t=mp.findroot(lambda z:mp_state(k,z,False)['Gt'],(tseed*mp.mpf('.99'),tseed*mp.mpf('1.01')),
        tol=mp.mpf(10)**(-dps+15),verify=True)
    R=mp_state(k,t,True,second)
    require(t>0 and abs(R['Gt'])<mp.mpf('1e-60') and R['Gtt']<0 and R['gap']>0,'Stationary local-maximum check failed.')
    require(mp_state(k,t*(1-mp.mpf('1e-4')),False)['Gt']>0 and mp_state(k,t*(1+mp.mpf('1e-4')),False)['Gt']<0,'Stationary sign bracket failed.')
    # Independent spectral exponential and divided-difference gain partial.
    M,B=mp_matrix(k);ev,W=mp.eig(M);Wi=W**-1
    Ts=W*mp.diag([mp.exp(z*t) for z in ev])*Wi;Z=Wi*B*W;F=mp.zeros(3)
    for i in range(3):
        for j in range(3):
            if i==j:f=t*mp.exp(ev[i]*t)
            else:
                require(abs(ev[i]-ev[j])>mp.mpf('1e-20'),'Spectral crosscheck near multiple roots.')
                f=mp.exp(ev[j]*t)*mp.expm1((ev[i]-ev[j])*t)/(ev[i]-ev[j])
            F[i,j]=f*Z[i,j]
    Ds=(R['u'].T*W*F*Wi*R['v'])[0]
    eT=mp.norm(Ts-R['T'])/max(1,mp.norm(R['T']))
    eD=abs(Ds-R['Gk'])/max(1,abs(R['Gk']))
    require(max(eT,eD)<mp.mpf('1e-55'),'Independent spectral comparison failed.')
    show=lambda z:mp.nstr(z,dps-10)
    row=dict(name=name,gain=show(k),time=show(t),peak=show(R['G']),alpha=show(max(mp.re(z) for z in ev)),
        gain_partial=show(R['Gk']),time_curvature=show(R['Gtt']),singular_gap=show(R['gap']),
        stationarity_residual=show(abs(R['Gt'])),spectral_propagator_error=show(eT),spectral_partial_error=show(eD),status=STATUS)
    if second:
        require(R['branch_second_derivative']>0,'Nonpositive local gain-branch curvature.')
        for n in ['Gkk','Gtk','branch_second_derivative']:row[n]=show(R[n])
    return row

def run(output:Path):
    output.parent.mkdir(parents=True,exist_ok=True)
    weak=float((-17+3*np.sqrt(217))/26);strong=float(963/320+9*np.sqrt(609)/64)
    # Direct time-envelope scans and local searches on the certified compact
    # retained intervals. Endpoints are explicit candidates, never omitted.
    double_runs=[];scanrows=[]
    for density in (1,2):
        for name,left in [('eta_half',weak),('eta_strong',strong)]:
            grid=np.linspace(left,10,65 if density==1 else 129)
            rows=[grid_peak(float(k),density,24.*density) for k in grid]
            for r in rows:scanrows.append(dict(case=name,**r))
            j=int(np.argmin([r['peak'] for r in rows]));lo=grid[max(0,j-1)];hi=grid[min(len(grid)-1,j+1)]
            opt=minimize_scalar(lambda k:grid_peak(float(k),density,24.*density)['peak'],bounds=(lo,hi),method='bounded',options={'xatol':1e-11})
            require(opt.success,'Gain search failed.')
            candidates=[rows[0],rows[-1],rows[j],grid_peak(float(opt.x),density,24.*density)]
            chosen=min(candidates,key=lambda r:r['peak'])
            double_runs.append(dict(case=name,left=left,right=10.,density=density,gain_samples=len(grid),
               selected=chosen,endpoint_selected=bool(chosen['gain']==left),
               local_optimizer=dict(gain=float(opt.x),peak=float(opt.fun)),
               sampled_curve_minimum_at_left=bool(j==0),
               scope='finite sampled search and local refinement; no global gain optimum certified'))
            print(f"double grid {density} {name}: k={chosen['gain']:.12g}, peak={chosen['peak']:.12g}",flush=True)
    # Detail samples resolve the shallow response difference.
    detail=[grid_peak(float(k)) for k in np.unique(np.r_[np.linspace(5,12,141),np.linspace(6.4,6.6,41),strong,6.44533,6.5,10,12])]
    # Spectrum samples use direct eigenvalues.
    spectral=[]
    for k in np.unique(np.r_[np.geomspace(.21,1000,241),np.linspace(5,40,141),weak,strong,6.44533,6.5,8,10,12,147/4]):
        eig=eigvals(matrix_float(float(k)));alpha=float(max(eig.real))
        polynomial=np.array([1.,k+5,7.5*k,13.5*k]);other=np.roots(polynomial)
        error=abs(alpha-float(max(other.real)))
        require(error<1e-9,'Spectral/cubic comparison failed.')
        spectral.append(dict(gain=float(k),alpha=alpha,cubic_alpha_error=error))
    runs=[]
    for dps in (80,110):
        mp.mp.dps=dps
        # Joint stationarity determines a candidate, not a global optimum.
        seed=double_runs[0]['selected'];seedt=mp.mpf(str(seed['time']))
        def joint(t,k):
            R=mp_state(k,t);return R['Gt'],R['Gk']
        optt,optk=mp.findroot(joint,(seedt,mp.mpf('6.44533')),tol=mp.mpf(10)**(-dps+15),maxsteps=25,verify=True)
        residual=max(abs(z) for z in joint(optt,optk));require(residual<mp.mpf('1e-60'),'Joint stationarity failed.')
        strongk=mp.mpf(963)/320+9*mp.sqrt(609)/64
        weakk=(-17+3*mp.sqrt(217))/26
        cases=[('weak_boundary',weakk),('low_gain',mp.mpf(2)),('peak_candidate',optk),
               ('strong_boundary',strongk),('strong_boundary_plus_1e-4',strongk+mp.mpf('0.0001')),
               ('incumbent',mp.mpf('6.5')),('exclusion_boundary',mp.mpf(10)),('reference_gain_12',mp.mpf(12))]
        rows=[high_precision_row(name,k,dps,optt if name=='peak_candidate' else None,name=='peak_candidate') for name,k in cases]
        by={r['name']:r for r in rows}
        require(mp.mpf(by['strong_boundary']['gain_partial'])>0,'Boundary right-branch partial is not positive.')
        require(mp.mpf(by['peak_candidate']['alpha'])>-mp.mpf(59)/20,'Peak candidate not excluded by stronger requirement.')
        cost=100*(mp.mpf(by['reference_gain_12']['peak'])/mp.mpf(by['peak_candidate']['peak'])-1)
        runs.append(dict(dps=dps,rows=rows,joint_stationarity_residual=mp.nstr(residual,30),
          relative_peak_cost_percent=mp.nstr(cost,dps-10),
          boundary_penalty_percent=mp.nstr(100*(mp.mpf(by['strong_boundary']['peak'])/mp.mpf(by['peak_candidate']['peak'])-1),dps-10)))
        print(f'{dps} digits: all {len(rows)} stationary and spectral checks pass',flush=True)
    comparisons=[]
    for low,high in zip(runs[0]['rows'],runs[1]['rows']):
        for key in ['gain','time','peak','alpha','gain_partial','time_curvature','singular_gap']:
            err=abs(mp.mpf(low[key])-mp.mpf(high[key]))/max(1,abs(mp.mpf(high[key])))
            require(err<mp.mpf('1e-55'),'Matched precision failed: '+key)
            comparisons.append(dict(name=high['name'],quantity=key,scaled_absolute_error=mp.nstr(err,25)))
    double_checks=[]
    hi={r['name']:r for r in runs[-1]['rows']}
    for name,row in hi.items():
        dense=grid_peak(float(row['gain']),2,48.)
        diff=abs(dense['peak']-float(row['peak']));require(diff<2e-12,'High-precision/time-grid mismatch.')
        double_checks.append(dict(name=name,height_difference=diff,grid=dense))
    require(all(x['endpoint_selected'] for x in double_runs if x['case']=='eta_strong'),'Constrained endpoint not selected.')
    for case in ['eta_half','eta_strong']:
        selected=[x['selected'] for x in double_runs if x['case']==case]
        require(abs(selected[0]['peak']-selected[1]['peak'])<2e-12,'Grid refinement mismatch.')
    obj=dict(status='passed',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      methods=dict(dps=[80,110],time_search='direct expm, mixed log/uniform grid, all sampled local maxima refined',
       gain_search='65/129 points on each retained interval, local refinement, both endpoints compared',
       high_precision='joint time/gain stationarity; positive conditional gain curvature at unconstrained candidate',
       gain_partial='6x6 block-exponential Frechet partial',hessian='9x9 block exponential plus simple-eigenvalue derivatives',
       independent_check='spectral exponential and divided differences',
       assurance='numerical witnesses/candidates only; q(T)<1 is floating-point tail evidence'),
      double_searches=double_runs,high_precision_runs=runs,precision_comparisons=comparisons,
      direct_grid_crosschecks=double_checks,stationary_evaluations=16,matched_precision_comparisons=len(comparisons),
      independent_spectral_comparisons=32,detail_peak_rows=len(detail),spectral_rows=len(spectral),
      constrained_gain_grid_rows=len(scanrows),new_global_peak_optimizer_certificate=False,
      new_all_time_upper_certificate=False)
    output.write_text(json.dumps(obj,indent=2)+'\n')
    for suffix,rows in [('witnesses',runs[-1]['rows']),('gain_grid',scanrows),('detail',detail),('spectrum',spectral)]:
        fields=list(dict.fromkeys(key for row in rows for key in row))
        with output.with_name('control_'+suffix+'.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    return obj
