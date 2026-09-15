"""Fixed physical matrices and declared numerical comparison inputs."""
from __future__ import annotations
import json
from pathlib import Path
import mpmath as mp
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
INPUTS=json.loads((ROOT/'configs/fixed_inputs.json').read_text())
def fixed_model(metric: str):
    if metric not in INPUTS['metric_factors']:
        raise ValueError('Unknown fixed performance metric')
    if mp.mp.dps > 180:
        raise ValueError('Declared decimal inputs support evaluation through 180 digits')
    T=mp.diag([mp.mpf(x) for x in INPUTS['metric_factors'][metric]])
    Aphysical=mp.matrix(INPUTS['physical_drift'])
    bphysical=mp.matrix(INPUTS['physical_input'])
    cphysical=mp.matrix(INPUTS['physical_measurement'])
    A0=T*Aphysical*T**-1
    b=T*bphysical;c=T**-1*cphysical
    given=INPUTS['prediction_inputs'][metric]
    out={k:mp.mpf(given[k]) for k in ('d','gamma','chi','a','h','j')}
    out.update(A0=A0,Q=-b*c.T/out['d'],T=T,P=T.T*T,
               u=mp.matrix(given['u']),v=mp.matrix(given['v']))
    return out

def matrix_float(metric: str, gain: float, shift: float=0.):
    if metric not in INPUTS['metric_factors'] or isinstance(gain,bool) or not np.isfinite(gain) or gain<0 or not np.isfinite(shift):
        raise ValueError('Invalid metric, gain or shift')
    from fractions import Fraction
    cast=lambda x:float(Fraction(x))
    T=np.diag([cast(x) for x in INPUTS['metric_factors'][metric]])
    A0=np.array([[cast(x) for x in row] for row in INPUTS['physical_drift']])
    b=np.array([cast(x) for x in INPUTS['physical_input']]);c=np.array([cast(x) for x in INPUTS['physical_measurement']])
    return T@(A0+gain*np.outer(b,c))@np.linalg.inv(T)+shift*np.eye(3)
