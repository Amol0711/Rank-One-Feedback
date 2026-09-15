"""Exact rational validation of certificate records and decay feasibility."""
from __future__ import annotations
from fractions import Fraction as F

def rational(value):
    if isinstance(value,bool) or not isinstance(value,(int,str,F)):
        raise TypeError('An exactly specified rational is required.')
    return F(value)

def shifted_cubic(gain,eta):
    k,e=rational(gain),rational(eta)
    if k<0 or e<0:raise ValueError('gain and margin must be nonnegative')
    return (F(1),5+k-3*e,3*e*e-10*e+(F(15,2)-2*e)*k,
            e*e*(5-e)+k*(e*e-F(15,2)*e+F(27,2)))

def strict_margin_record(gain,eta):
    p=shifted_cubic(gain,eta);q=p[1]*p[2]-p[3]
    vals={'quadratic_coefficient':p[1],'linear_coefficient':p[2],
          'constant_coefficient':p[3],'routh_determinant':q}
    return {'gain':str(gain),'eta':str(eta),
            'shifted_characteristic_coefficients':[str(z) for z in p],
            'positive_conditions':{key:str(z) for key,z in vals.items()},
            'strictly_feasible':all(z>0 for z in vals.values()),
            'meaning':'spectral abscissa strictly below negative eta'}

def verify_partition(rows,expected_cells):
    if isinstance(expected_cells,bool) or not isinstance(expected_cells,int) or expected_cells<1:
        raise ValueError('positive cell count required')
    if len(rows)!=expected_cells:raise ValueError('incomplete cell set')
    for j,row in enumerate(rows):
        if row['cell']!=j or F(row['x_left'])!=F(j,expected_cells) or F(row['x_right'])!=F(j+1,expected_cells):
            raise ValueError('cell gap, overlap, wrong endpoint, or wrong ordering')
    return True

def endpoint_fraction(endpoint):
    sign,m,e,bc=map(int,endpoint)
    if sign not in (0,1) or m<0 or bc<0:raise ValueError('nonfinite or malformed endpoint')
    q=F((-1 if sign else 1)*m)
    return q*2**e if e>=0 else q/F(2**(-e))

def verify_serialized_interval(record):
    lo,hi=map(endpoint_fraction,record['binary_endpoints'])
    if not F(record['lower'])<=lo<=hi<=F(record['upper']):
        raise ValueError('decimal serialization is not outward')
    return lo,hi

def verify_record_structure(record,case):
    u=record['upper_certificate'];l=record['exclusion_certificate']
    for key in ('metric',):
        if u[key]!=case[key] or l[key]!=case[key]:raise ValueError('metric mismatch')
    for key in ('gain','step'):
        if F(u[key])!=F(case[key]):raise ValueError('upper input mismatch')
    for key in ('cutoff','tau'):
        if F(l[key])!=F(case[key]):raise ValueError('lower input mismatch')
    if l['u_rational']!=case['u_rational'] or l['v_rational']!=case['v_rational']:
        raise ValueError('witness mismatch')
    if l['cells']!=case['cells']:raise ValueError('cell count mismatch')
    verify_partition(l['cell_records'],case['cells'])
    if u.get('includes_time_zero') is not True:raise ValueError('time zero is not explicitly covered')
    if u['steps']<1 or len(u['samples'])!=u['steps']:
        raise ValueError('missing positive-time samples')
    for j,row in enumerate(u['samples'],1):
        if row['step']!=j:raise ValueError('sample gap or wrong ordering')
        verify_serialized_interval(row['norm_upper'])
    for label in ('upper','grid_upper','mu_upper','stop_time','stop_norm_upper'):
        verify_serialized_interval(u[label])
    left,right=verify_serialized_interval(u['stop_time'])
    if not left<=F(u['step'])*u['steps']<=right:raise ValueError('wrong stop time')
    if F(u['stop_norm_upper']['upper'])>=1:raise ValueError('no strict tail stop')
    sample_max=max(F(1),*(endpoint_fraction(x['norm_upper']['binary_endpoints'][1]) for x in u['samples']))
    if endpoint_fraction(u['grid_upper']['binary_endpoints'][1])<sample_max:raise ValueError('grid maximum does not cover every sample')
    for row in l['cell_records']:verify_serialized_interval(row['scalar_witness'])
    verify_serialized_interval(l['lower'])
    if endpoint_fraction(l['lower']['binary_endpoints'][0])>min(endpoint_fraction(x['scalar_witness']['binary_endpoints'][0]) for x in l['cell_records']):raise ValueError('aggregate lower bound does not cover every cell')
    if not (F(u['upper']['upper'])<F(case['display_upper'])<F(case['display_lower'])<F(l['lower']['lower'])):
        raise ValueError('no strict incumbent/exclusion separation')
    if not strict_margin_record(case['gain'],case['eta'])['strictly_feasible']:
        raise ValueError('incumbent fails strict decay margin')
    return True
