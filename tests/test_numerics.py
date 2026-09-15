"""Regression and adverse-input tests for numerical calculations and records."""
from __future__ import annotations
import copy,json,hashlib
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
import numpy as np
import pytest
from numerics.fixed_model import fixed_model,matrix_float,INPUTS,ROOT
from numerics import interval,integer_replay,stationary,directions
from numerics.certificate_checks import verify_record_structure,verify_partition,verify_serialized_interval,strict_margin_record
from numerics.workflows import selected_rows,check_reference,numeric_equal
from reproduce import destination
from verify_release import verify

@pytest.mark.parametrize('metric',['raw','normalized'])
@pytest.mark.parametrize('gain',[0,0.354,1.046,6.5,10,1000])
def test_fixed_metric_similarity(metric,gain):
    with mp.workdps(80):
        m=fixed_model(metric)
        expected=np.array((m['A0']-mp.mpf(str(gain))*m['Q']).tolist(),float)
        assert np.max(abs(matrix_float(metric,gain)-expected))<1e-12
        assert abs(float(m['d'])-1)<1e-14

@pytest.mark.parametrize('metric',['raw','normalized'])
def test_unit_linked_inputs(metric):
    with mp.workdps(110):
        m=fixed_model(metric)
        assert abs(mp.norm(m['u'])-1)<mp.mpf('1e-105')
        assert abs(mp.norm(m['v'])-1)<mp.mpf('1e-105')
        assert m['gamma']>1 and m['chi']>0 and m['a']>0

@pytest.mark.parametrize('value',[-1,float('nan'),float('inf'),True])
def test_bad_gain(value):
    with pytest.raises(ValueError):matrix_float('normalized',value)

def test_bad_metric():
    with pytest.raises(ValueError):fixed_model('other')

def test_precision_input_limit():
    with mp.workdps(181),pytest.raises(ValueError):fixed_model('normalized')

@pytest.mark.parametrize('metric',['raw','normalized'])
def test_stationary_fixed_case(metric):
    with mp.workdps(80):
        r=stationary.evaluate(metric,50,80)
        ref=json.loads((ROOT/'reference/stationary.json').read_text())['runs'][0]['rows']
        old=next(x for x in ref if x['metric']==metric and x['gain']==50)
        for k in ('peak','peak_time','gain_derivative','R_G','R_t','R_D'):numeric_equal(old[k],r[k])
        assert mp.mpf(r['curvature'])<0 and mp.mpf(r['singular_gap'])>0

@pytest.mark.parametrize('delta',[0,1])
def test_shift_and_linked_direction(delta):
    r=directions.compute_case(50,delta,80)
    assert r['numerical_checks_passed']
    with mp.workdps(100):
        for k in ('physical_initial_P_norm_residual','physical_response_P_norm_residual','direction_transport_residual','drift_shift_identity_residual'):
            assert mp.mpf(r['checks'][k])<mp.mpf('1e-55')

@pytest.mark.parametrize('value',[0.5,True,float('nan')])
def test_interval_reject_binary_inputs(value):
    with pytest.raises(TypeError):interval.exact(value)

@pytest.mark.parametrize('digits',[0,10,29,True,31.5])
def test_interval_precision_validation(digits):
    with pytest.raises(ValueError):interval.precision(digits)

@pytest.mark.parametrize('value',['-3/7','1/3','0','19/13'])
def test_outward_serialization(value):
    interval.precision(50);r=interval.display_bounds(interval.exact(value),7)
    lo,hi=verify_serialized_interval(r)
    assert F(r['lower'])<=lo<=F(value)<=hi<=F(r['upper'])

@pytest.mark.parametrize('cells',[1,3,128])
def test_exact_partition(cells):
    rows=[dict(cell=j,x_left=f'{j}/{cells}',x_right=f'{j+1}/{cells}') for j in range(cells)]
    assert verify_partition(rows,cells)
    rows[-1]['x_right']='2'
    with pytest.raises(ValueError):verify_partition(rows,cells)

@pytest.mark.parametrize('metric,eta,gain',[('raw','0.1','0.354'),('raw','0.5','1.046'),('normalized','0.5','6.5'),('normalized','59/20','6.5')])
def test_fixed_incumbent_feasibility(metric,eta,gain):
    assert strict_margin_record(gain,eta)['strictly_feasible']

@pytest.mark.parametrize('digits',[50,80])
@pytest.mark.parametrize('case_index',[0,1,2])
def test_certificate_complete(digits,case_index):
    case=json.loads((ROOT/'configs/certificate_cases.json').read_text())[case_index]
    r=json.loads((ROOT/f'reference/interval_{digits}/{case["name"]}.json').read_text())
    assert verify_record_structure(r,case)

@pytest.mark.parametrize('mutation',['sample','cell','witness','metric','stop','upper','precision_serialization','time_zero','grid'])
def test_certificate_mutations(mutation):
    case=json.loads((ROOT/'configs/certificate_cases.json').read_text())[2]
    r=json.loads((ROOT/f'reference/interval_50/{case["name"]}.json').read_text())
    if mutation=='time_zero':r['upper_certificate']['includes_time_zero']=False
    if mutation=='grid':r['upper_certificate']['grid_upper']=interval.display_bounds(interval.exact('1.0'))
    if mutation=='sample':r['upper_certificate']['samples'].pop(0)
    if mutation=='cell':r['exclusion_certificate']['cell_records'].pop()
    if mutation=='witness':r['exclusion_certificate']['u_rational'][0]='0'
    if mutation=='metric':r['upper_certificate']['metric']='raw'
    if mutation=='stop':r['upper_certificate']['stop_norm_upper']['upper']='1.0'
    if mutation=='upper':r['upper_certificate']['upper']['upper']='1.5'
    if mutation=='precision_serialization':r['exclusion_certificate']['cell_records'][0]['scalar_witness']['lower']='1000'
    with pytest.raises(ValueError):verify_record_structure(r,case)

@pytest.mark.parametrize('a,b',[('1/3','-2/7'),('-3/4','5/11'),('1','2'),('0','1/3')])
def test_integer_enclosing_operations(a,b):
    x,y=integer_replay.DI.scalar(a),integer_replay.DI.scalar(b)
    for z,q in [(x+y,F(a)+F(b)),(x-y,F(a)-F(b)),(x*y,F(a)*F(b)),(x/y,F(a)/F(b))]:
        lo,hi=z.endpoints();assert lo<=q<=hi

def test_integer_reject_zero_denominator():
    with pytest.raises(ZeroDivisionError):integer_replay.DI.scalar(1)/integer_replay.DI.range(-1,1)

def test_selection_is_single_source():
    spec=json.loads((ROOT/'configs/diagnostic_selection.json').read_text());rows=selected_rows()
    assert len(rows)==len(spec['rows'])==8
    assert spec['schema_version']==2 and spec['display_columns']==spec['columns'][1:]
    assert all(r['metric']=='normalized' for r in rows)
    assert [int(r['gain']) for r in rows]==[10,25,50,100,140,150,1000,2**35]
    assert [(r['metric'],int(r['gain'])) for r in rows]==[(r['metric'],r['gain']) for r in spec['rows']]
    assert int(rows[-1]['gain'])==2**35 and rows[-1]['correction_error_pct'] is None
    assert all(int(r['gain'])!=400 for r in rows)
    full=list(__import__('csv').DictReader((ROOT/'data/asymptotic_comparison.csv').open()))
    assert len([r for r in full if r['metric']=='normalized' and r['gain']=='400'])==1

def test_stored_reference_checks():assert check_reference()['passed']

def test_package_hashes():assert verify()['passed']

@pytest.mark.parametrize('change',['modify','extra','missing','nested_manifest','symlink'])
def test_manifest_rejects_mutations(tmp_path,change):
    p=tmp_path/'copy';p.mkdir();f=p/'value.txt';f.write_text('1')
    m={'value.txt':{'bytes':1,'sha256':hashlib.sha256(b'1').hexdigest()}}
    (p/'MANIFEST.json').write_text(json.dumps({'files':m}))
    if change=='modify':f.write_text('2')
    if change=='extra':(p/'extra.txt').write_text('1')
    if change=='missing':f.unlink()
    if change=='nested_manifest':
        (p/'extra').mkdir();(p/'extra/MANIFEST.json').write_text('{}')
    if change=='symlink':(p/'link').symlink_to(f)
    with pytest.raises(ValueError):verify(p)

@pytest.mark.parametrize('kind',['existing','internal','ancestor','symlink'])
def test_output_rejects_unsafe_paths(tmp_path,kind):
    if kind=='existing':path=tmp_path
    elif kind=='internal':path=ROOT/'unused_output'
    elif kind=='ancestor':path=ROOT.parent
    else:
        link=tmp_path/'link';link.symlink_to(tmp_path,target_is_directory=True);path=link/'out'
    with pytest.raises(ValueError):destination(path)

def test_new_external_output(tmp_path):assert destination(tmp_path/'new')==tmp_path/'new'


def test_local_checkout_metadata_is_ignored(tmp_path):
    p=tmp_path/'root';p.mkdir();(p/'MANIFEST.json').write_text(json.dumps({'files':{}}))
    (p/'.git').mkdir();(p/'.git/HEAD').write_text('local checkout bookkeeping')
    assert verify(p)['passed']


@pytest.mark.parametrize('mutation',['duplicate','raw','missing','columns','omission','source','version'])
def test_selection_contract_rejects_mutation(monkeypatch,mutation):
    from numerics import workflows
    original=workflows.load
    spec=copy.deepcopy(original(ROOT/'configs/diagnostic_selection.json'))
    if mutation=='duplicate':spec['rows'][1]=copy.deepcopy(spec['rows'][0])
    if mutation=='raw':spec['rows'][0]['metric']='raw'
    if mutation=='missing':spec['rows'].pop()
    if mutation=='columns':spec['display_columns'].pop()
    if mutation=='omission':spec['omission_token']='0'
    if mutation=='source':spec['rows'][0]['source']='unknown'
    if mutation=='version':spec['schema_version']=0
    monkeypatch.setattr(workflows,'load',lambda p: spec if Path(p).name=='diagnostic_selection.json' else original(p))
    with pytest.raises(ValueError):workflows.selected_rows()
