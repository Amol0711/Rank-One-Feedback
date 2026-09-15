"""Numerical workflow dispatch and comparisons with retained reference values."""
from __future__ import annotations
import csv,io,json,hashlib,time
from pathlib import Path
from decimal import Decimal,localcontext
from fractions import Fraction
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from .fixed_model import ROOT,INPUTS,matrix_float
from . import stationary,directions,finite_gain,interval,integer_replay
from .certificate_checks import verify_record_structure,strict_margin_record,endpoint_fraction

MAIN_QUANTITIES=('peak','peak_time','gain_derivative','predicted_peak','predicted_time','predicted_gain_derivative','R_G','R_t','R_D')

def require(ok,message):
    if not ok:raise ValueError(message)

def load(path):return json.loads(Path(path).read_text())
def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def csv_read(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def csv_write(path,rows,fields=None):
    require(bool(rows),'Empty numerical row set')
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def numeric_equal(a,b,tolerance='1e-55',floor='1e-100'):
    with mp.workdps(160):
        x,y=mp.mpf(a),mp.mpf(b)
        require(mp.isfinite(x) and mp.isfinite(y),'Nonfinite reference comparison')
        require(abs(x-y)<=mp.mpf(tolerance)*max(abs(x),abs(y),mp.mpf(floor)),f'Numerical reference mismatch for {a} and {b}')

def compare_rows(left,right,keys,identity=('metric','gain','dps'),tolerance='1e-55'):
    require(len(left)==len(right),'Reference row coverage differs')
    checks=0
    for a,b in zip(left,right):
        for k in identity:require(a[k]==b[k],'Case identity differs '+k)
        for k in keys:
            if a[k] is None or isinstance(a[k],bool):require(a[k]==b[k],'Omission/flag mismatch')
            else:numeric_equal(a[k],b[k],tolerance)
            checks+=1
    return checks

def selected_rows():
    spec=load(ROOT/'configs/diagnostic_selection.json')
    require(spec['schema_version']==2 and spec['selection_id']=='diagnostic_selection_v2',
            'Unsupported selected-view version')
    require(spec['columns']==['metric','gain','height_error_pct','correction_error_pct',
            'time_error_pct','sensitivity_error_pct','R_G'] and
            spec['display_columns']==spec['columns'][1:], 'Selected-view columns differ')
    require(spec['omission_token']=='' and spec['correction_floor_relative_to_gamma']=='1e-8',
            'Selected-view omission convention differs')
    identities=[(r['metric'],r['gain']) for r in spec['rows']]
    require(len(identities)==8 and len(set(identities))==8 and
            all(m=='normalized' and type(g) is int and g>0 for m,g in identities) and
            [g for m,g in identities]==sorted(g for m,g in identities),
            'Selected-view identities differ')
    all_rows=csv_read(ROOT/spec['sources']['grid'])
    special=load(ROOT/spec['sources']['large_gain'])['runs'][-1]['row']
    out=[]
    for item in spec['rows']:
        if item['source']=='grid':
            matches=[r for r in all_rows if r['metric']==item['metric'] and Decimal(r['gain'])==item['gain']]
            require(len(matches)==1,'Selected numerical row is missing or duplicated');r=matches[0]
        elif item['source']=='large_gain':
            r=special;require(r['gain']==item['gain'] and r['metric']==item['metric'],'Special case identity differs')
        else:raise ValueError('Unknown numerical selection source')
        out.append({k:r[k] for k in spec['columns']})
    return out

def check_reference():
    s=load(ROOT/'reference/stationary.json')
    require(s['gains']==list(stationary.GAINS) and s['metrics']==list(stationary.METRICS),'Stationary input grid differs')
    require(len(s['runs'])==2,'Both precision runs are required')
    checks=0
    for run in s['runs']:
        require(len(run['rows'])==44,'Incomplete fixed-metric grid')
        for r in run['rows']:
            with mp.workdps(140):
                vals=[mp.mpf(r[k]) for k in ('peak','peak_time','gain_derivative','predicted_peak','predicted_time','predicted_gain_derivative','gamma','r','L','d')]
                expected=stationary.discrepancy_values(*vals)
                for k,v in expected.items():
                    if v is None or isinstance(v,bool):require(r[k]==v,'Omission/flag mismatch '+k)
                    else:numeric_equal(r[k],v,'1e-45');checks+=1
                require(mp.mpf(r['curvature'])<0 and mp.mpf(r['singular_gap'])>0,'Stationary curvature/gap invalid')
                require(mp.mpf(r['stationarity_residual'])<mp.mpf('1e-60'),'Stationary residual invalid')
    require(len(csv_read(ROOT/'data/asymptotic_comparison.csv'))==44,'CSV coverage mismatch')
    require(len(csv_read(ROOT/'data/peak_diagnostics.csv'))==22,'Normalized CSV coverage mismatch')
    certs=[]
    for digits in (50,80):
        for case in load(ROOT/'configs/certificate_cases.json'):
            r=load(ROOT/f'reference/interval_{digits}/{case["name"]}.json')
            verify_record_structure(r,case)
            u=r['upper_certificate'];lo=r['exclusion_certificate']
            require(u['includes_time_zero'] is True,'Time zero omitted')
            require(endpoint_fraction(u['grid_upper']['binary_endpoints'][1])>=max(endpoint_fraction(x['norm_upper']['binary_endpoints'][1]) for x in u['samples']),'Grid upper smaller than a sample')
            require(endpoint_fraction(lo['lower']['binary_endpoints'][0])<=min(endpoint_fraction(x['scalar_witness']['binary_endpoints'][0]) for x in lo['cell_records']),'Aggregate lower exceeds a cell bound')
            certs.append({'case':case['name'],'digits':digits,'samples':u['steps'],'cells':lo['cells']})
    # The tighter spectral condition reuses the same normalized peak calculations.
    require(strict_margin_record('6.5','59/20')['strictly_feasible'],'Stronger incumbent is infeasible')
    sel=selected_rows();require(len(sel)==len(load(ROOT/'configs/diagnostic_selection.json')['rows']),'Selection coverage differs')
    return dict(passed=True,stored_numeric_comparisons=checks,certificate_records=certs,selected_rows=len(sel),kind='stored-record arithmetic and schema checks; no matrix enclosure or search recomputed')

def export(out):
    for src in sorted((ROOT/'data').glob('*.csv')):(out/src.name).write_bytes(src.read_bytes())
    csv_write(out/'selected_diagnostics.csv',selected_rows())
    return dict(passed=True,csv_files=16,kind='retained numerical data export; no solver invoked')

def responses(out):
    total=0;maximum_error=0.
    for name,spec in load(ROOT/'configs/response_grids.json').items():
        gain=spec['gain'];t=np.array([float(z) for z in spec['times']]);A=matrix_float('normalized',gain)
        v=np.linalg.svd(expm(t[:,None,None]*A),compute_uv=False)[:,0]
        ref=csv_read(ROOT/'data'/f'{name}.csv');require(len(ref)==len(t),'Time grid changed')
        rows=[]
        for i,(ti,vi) in enumerate(zip(t,v)):
            er=abs(float(vi)-float(ref[i]['norm_witness']));maximum_error=max(maximum_error,er)
            require(er<=2e-10*max(1,abs(float(vi))),'Response sample mismatch')
            with localcontext() as ctx:
                ctx.prec=80;fast=str(Decimal(spec['times'][i])*gain)
            rows.append(dict(time=spec['times'][i],fast_time=fast,norm_witness=repr(float(vi))))
        csv_write(out/f'{name}.csv',rows);total+=len(rows)
    return dict(passed=True,response_samples=total,max_absolute_reference_difference=maximum_error,kind='fresh matrix-exponential samples at fixed input times')

def stationary_grid(out):
    r=stationary.run(out/'stationary.json');old=load(ROOT/'reference/stationary.json');checks=0
    for a,b in zip(old['runs'],r['runs']):
        checks+=compare_rows(a['rows'],b['rows'],MAIN_QUANTITIES+('correction_error_pct','correction_defined','sensitivity_error_pct','log_parameter_below_one','below_supplied_comparison_range'))
    return dict(passed=True,evaluations=r['stationary_evaluations'],reference_comparisons=checks,kind='fresh stationary search at 80 and 110 digits; not interval-certified extrema')

def direction_grid(out):
    r=directions.run(out);old=load(ROOT/'reference/directions.json');checks=0
    require(len(r['rows'])==len(old['rows'])==11,'Direction case coverage differs')
    for a,b in zip(old['rows'],r['rows']):
        require((a['gain'],a['delta'])==(b['gain'],b['delta']),'Direction case differs')
        for prec in ('at_80_digits','at_110_digits'):
            x,y=a[prec],b[prec]
            for k in x['values']:numeric_equal(x['values'][k],y['values'][k]);checks+=1
            for k in x['vectors']:
                for u,v in zip(x['vectors'][k],y['vectors'][k]):numeric_equal(u,v);checks+=1
            require(x['in_supplied_comparison_range']==y['in_supplied_comparison_range'],'Range flag differs')
    return dict(passed=True,cases=11,reference_comparisons=checks,kind='fresh linked-direction and independently shifted stationary calculations')

def comparison(out):
    inp=INPUTS['comparison_inputs'];gain=inp['rate'];old=load(ROOT/'reference/comparison_residual.json');runs=[];checks=0
    for digits in (80,110):
        with mp.workdps(digits):
            row=stationary.evaluate('normalized',gain,digits)
            rr,L,d=[mp.mpf(row[k]) for k in ('r','L','d')]
            scales={'height':inp['height_multiplier']*L**2/rr**2,'time':inp['time_multiplier']*L/rr**2,'sensitivity':inp['sensitivity_multiplier']*d*L**2/rr**3}
            comps={}
            for name,bound in scales.items():
                residual=abs(mp.mpf(row[name+'_residual']));require(0<residual<bound,'Residual comparison failed')
                comps[name]={k:mp.nstr(v,digits-10) for k,v in dict(absolute_residual=residual,applicable_bound=bound,bound_to_observed_ratio=bound/residual).items()}
            runs.append(dict(dps=digits,row=row,bound_comparisons=comps))
    for a,b in zip(old['runs'],runs):
        checks+=compare_rows([a['row']],[b['row']],MAIN_QUANTITIES+('correction_error_pct','correction_defined'))
        for n in a['bound_comparisons']:
            for k,v in a['bound_comparisons'][n].items():numeric_equal(v,b['bound_comparisons'][n][k]);checks+=1
    save(out/'comparison_residual.json',dict(runs=runs,supplied_inputs=inp,assurance='Stationary numerical residuals compared with fixed supplied scales; no interval enclosure of residuals.'))
    return dict(passed=True,evaluations=2,reference_comparisons=checks,correction_error_omitted=True,comparison_inputs_recomputed=False)

def finite_search(out):
    r=finite_gain.run(out/'finite_gain.json');old=load(ROOT/'reference/finite_gain.json');checks=0
    for a,b in zip(old['high_precision_runs'],r['high_precision_runs']):
        require(len(a['rows'])==len(b['rows']),'Finite-gain row count changed')
        for x,y in zip(a['rows'],b['rows']):
            require(x['name']==y['name'],'Finite-gain row identity changed')
            for key in ('gain','time','peak','alpha','gain_partial','time_curvature','singular_gap'):
                with mp.workdps(150):
                    if key=='gain_partial' and abs(mp.mpf(x[key]))<mp.mpf('1e-60'):
                        require(abs(mp.mpf(x[key])-mp.mpf(y[key]))<mp.mpf('1e-60'),'Near-zero gain partial differs')
                    else:numeric_equal(x[key],y[key],'1e-50')
                checks+=1
    return dict(passed=True,stationary_evaluations=r['stationary_evaluations'],gain_grid_rows=r['constrained_gain_grid_rows'],detail_rows=r['detail_peak_rows'],spectral_rows=r['spectral_rows'],reference_comparisons=checks,global_optimizer_certified=False)

def intervals(out):
    summary=[]
    for digits in (50,80):
        interval.precision(digits);dest=out/f'interval_{digits}';dest.mkdir()
        for case in load(ROOT/'configs/certificate_cases.json'):
            up=interval.fixed_gain_upper(case['gain'],case['step'],case['metric'])
            lo=interval.high_gain_witness(case['cutoff'],case['tau'],case['metric'],case['cells'],(case['u_rational'],case['v_rational']))
            r=dict(name=case['name'],upper_certificate=up,exclusion_certificate=lo,display_upper=case['display_upper'],display_lower=case['display_lower'],strict_margin=strict_margin_record(case['gain'],case['eta']))
            verify_record_structure(r,case);save(dest/f'{case["name"]}.json',r)
            # Numerical interval endpoints may vary with the proposal backend.
            # Recover the same strict displayed inequalities rather than demand identical endpoints.
            summary.append(dict(case=case['name'],digits=digits,upper=up['upper']['upper'],lower=lo['lower']['lower'],samples=up['steps'],cells=lo['cells']))
            print('interval',digits,case['name'],'PASS',flush=True)
    return dict(passed=True,cases=summary,kind='fresh directed interval matrix exponentials, all positive-time samples and every reciprocal-gain cell')

def replay(out):
    rows=[]
    for case in load(ROOT/'configs/certificate_cases.json'):
        r=integer_replay.replay_case(load(ROOT/f'reference/interval_50/{case["name"]}.json'),case)
        save(out/f'{case["name"]}.json',r);rows.append({k:r[k] for k in ('name','upper','lower','positive_time_samples')})
        print('integer replay',case['name'],'PASS',flush=True)
    return dict(passed=True,cases=rows,kind='independent integer-grid re-evaluation using retained rational proposals')

WORKFLOWS={'export':export,'responses':responses,'stationary':stationary_grid,'directions':direction_grid,'comparison':comparison,'finite-gain':finite_search,'intervals':intervals,'integer-replay':replay}
