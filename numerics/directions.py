"""Linked singular vectors and shifted stationary responses at fixed gains."""
from __future__ import annotations
import csv,json,math
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import expm,svdvals
from scipy.optimize import minimize_scalar
from .fixed_model import fixed_model,INPUTS
DIRECTION_GAINS=(50,100,1000,10000,2**35)
WEIGHTED_GAINS=(10,50,100,1000,100000,2**38)
PRECISIONS=(80,110)
R=INPUTS['comparison_inputs']['rate']
R_DELTA=INPUTS['comparison_inputs']['shifted_rate']
def model(gain: int, delta: int=0):
    if isinstance(gain,bool) or gain<=0 or int(gain)!=gain or delta not in (0,1):
        raise ValueError('Positive integer gain and shift 0 or 1 required')
    d=fixed_model('normalized');k=mp.mpf(gain)
    A=d['A0']-k*d['Q']+delta*mp.eye(3)
    Aph=d['T']**-1*A*d['T']
    return dict(k=k,A=A,Aph=Aph,TP=d['T'],P=d['P'],gamma=d['gamma'],chi=d['chi'],ap=d['a'],u=d['u'],v=d['v'])

def linked_pair(T, uref, vref, right_gram=False):
    """The companion vector is transported, never signed independently."""
    values, vectors = mp.eigsy(T.T*T if right_gram else T*T.T)
    s = mp.sqrt(values[2])
    if right_gram:
        v = vectors[:, 2]
        u = T*v/s
    else:
        u = vectors[:, 2]
        v = T.T*u/s
    overlap = (uref.T*u + vref.T*v)[0]
    if overlap == 0:
        raise ArithmeticError('Cannot choose the positive-overlap linked sign.')
    if overlap < 0:
        u, v = -u, -v
    return s, u, v, values, vectors

def pair_checks(T, s, u, v, uref, vref, TP, P, Aph, t):
    x, y = TP**-1*v, TP**-1*u
    xref, yref = TP**-1*vref, TP**-1*uref
    pnorm = lambda z: mp.sqrt((z.T*P*z)[0])
    return dict(
        u_norm_residual=abs(mp.norm(u)-1),
        v_norm_residual=abs(mp.norm(v)-1),
        linked_forward_residual=mp.norm(T*v-s*u),
        linked_adjoint_residual=mp.norm(T.T*u-s*v),
        linked_overlap=(uref.T*u+vref.T*v)[0],
        physical_initial_P_norm_residual=abs(pnorm(x)-1),
        physical_response_P_norm_residual=abs(pnorm(y)-1),
        physical_map_residual=pnorm(mp.expm(Aph*t)*x-s*y),
        transported_direction_error=pnorm(y-yref)+pnorm(x-xref),
        direction_transport_residual=abs(pnorm(y-yref)+pnorm(x-xref)-mp.norm(u-uref)-mp.norm(v-vref)),
    )

def coarse_search(gain: int, delta: int, density: int = 1):
    """Floating-point search used for starting values only, not a bound."""
    if gain >= R:
        raise ValueError('Use the high-precision fast-time seed at the supplied large-gain comparison points.')
    A0 = np.array([[0., 5., 0.], [0., 0., 0.], [.5, 0., -5.]])
    A = A0 + gain*np.outer([1., .6, 0.], [-1., 0., 1.]) + delta*np.eye(3)
    alpha = float(np.linalg.eigvals(A).real.max())
    if alpha >= 0:
        raise ArithmeticError('Non-Hurwitz matrix in the sampled case.')
    horizon = max(8., 24./(-alpha))*density
    times = np.unique(np.r_[0., np.geomspace(1e-7/gain, horizon, 700*density),
                            np.linspace(0, horizon, 350*density)])
    heights = np.array([svdvals(expm(A*t))[0] for t in times])
    ids = np.where((heights[1:-1] > heights[:-2]) &
                   (heights[1:-1] >= heights[2:]) & (heights[1:-1] > 1.00000001))[0]+1
    humps = []
    for i in ids:
        result = minimize_scalar(lambda tau: -svdvals(expm(A*(tau/gain)))[0],
                                 bounds=(times[i-1]*gain, times[i+1]*gain), method='bounded',
                                 options={'xatol': 1e-12})
        if not result.success:
            raise ArithmeticError('Seed search failed.')
        humps.append({'fast_time': float(result.x), 'height': float(-result.fun)})
    if not humps:
        raise ArithmeticError('No above-one numerical response hump found.')
    best = max(humps, key=lambda p: p['height'])
    return dict(gain=gain, delta=delta, density=density, horizon=horizon,
                sample_count=len(times), terminal_norm=float(heights[-1]),
                spectral_abscissa=alpha, humps=humps, selected=best,
                assurance='ordinary floating-point seed search; not an all-time enclosure')

def compute_case(gain: int, delta: int, dps: int, seed=None):
    with mp.workdps(dps):
        m = model(gain, delta)
        A, k = m['A'], m['k']
        B = A/k
        L = mp.log(m['chi']*k/m['ap'])
        Ld = mp.log(m['chi']*k/(m['ap']-delta*m['gamma']))

        def state(tau):
            T = mp.expm(B*tau)
            return T, linked_pair(T, m['u'], m['v'])

        def g(tau):
            _, (_, u, _, _, _) = state(tau)
            return (u.T*B*u)[0]  # derivative of log(sigma) in fast time

        start = mp.mpf(str(seed)) if seed is not None else Ld
        tau = mp.findroot(g, (start*mp.mpf('.99'), start*mp.mpf('1.01')),
                          tol=mp.mpf(10)**(-(dps-18)), maxsteps=80)
        if not tau > 0:
            raise ArithmeticError('Nonpositive stationary time.')
        t = tau/k
        T, (s, u, v, lam, U) = state(tau)
        err = mp.norm(u-m['u'])+mp.norm(v-m['v'])
        checks = pair_checks(T, s, u, v, m['u'], m['v'], m['TP'], m['P'], m['Aph'], t)
        H = T*T.T
        H1 = B*H+H*B.T
        H2 = B*B*H+2*B*H*B.T+H*B.T*B.T
        lam1 = (u.T*H1*u)[0]
        lam2 = (u.T*H2*u)[0]+2*sum((U[:,j].T*H1*u)[0]**2/(lam[2]-lam[j]) for j in (0,1))
        sigma2 = lam2/(2*s)-lam1**2/(4*s**3)
        stationarity = (u.T*B*u)[0]
        step = mp.mpf('1e-6')
        before, after = g(tau*(1-step)), g(tau*(1+step))
        evals, V = mp.eig(A, left=False, right=True)
        Tspectral = V*mp.diag([mp.exp(z*t) for z in evals])*V**-1
        sright, uright, vright, _, _ = linked_pair(T, m['u'], m['v'], right_gram=True)
        Tphysical = m['TP']*mp.expm(m['Aph']*t)*m['TP']**-1
        Tunshifted = mp.expm((A-delta*mp.eye(3))*t)
        checks.update(
            stationarity_residual_fast=abs(stationarity),
            stationarity_residual_time=abs(stationarity*k),
            stationarity_formula_residual=abs(lam1/(2*s*s)-stationarity),
            spectral_exponential_residual=mp.norm(Tspectral-T)/max(1,mp.norm(T)),
            right_gram_height_residual=abs(sright-s),
            right_gram_pair_residual=mp.norm(uright-u)+mp.norm(vright-v),
            physical_similarity_residual=mp.norm(Tphysical-T)/max(1,mp.norm(T)),
            drift_shift_identity_residual=mp.norm(T-mp.exp(delta*t)*Tunshifted)/max(1,mp.norm(T)),
            stationary_curvature_fast=sigma2,
            derivative_before_fast=before,
            derivative_after_fast=after,
            top_singular_gap=s-mp.sqrt(lam[1]),
        )
        threshold = mp.mpf(10)**(-(dps-25))
        residual_fields = [key for key in checks if 'residual' in key]
        passed = (all(checks[key] < threshold for key in residual_fields)
                  and checks['linked_overlap'] > 0 and sigma2 < 0 and before > 0 and after < 0
                  and checks['top_singular_gap'] > 0 and 1 < s < m['gamma'])
        if not passed:
            raise ArithmeticError(f'Failed checks at k={gain}, delta={delta}, dps={dps}: {checks}')
        values = dict(stationary_time=t, fast_time=tau, stationary_height=s, gamma=m['gamma'],
                      gap_below_gamma=m['gamma']-s, logarithmic_time=Ld)
        if delta == 0:
            values.update(direction_error=err, comparison_112L_over_r=INPUTS['comparison_inputs']['direction_multiplier']*L/k,
                          comparison_to_error_ratio=(INPUTS['comparison_inputs']['direction_multiplier']*L/k)/err,
                          scaled_direction_error=err*k/L)
        encode = lambda x: mp.nstr(x, dps-8)
        vectors = dict(u=[encode(x) for x in u], v=[encode(x) for x in v],
                       physical_initial=[encode(x) for x in m['TP']**-1*v],
                       physical_response=[encode(x) for x in m['TP']**-1*u])
        return dict(gain=gain, delta=delta, decimal_precision=dps,
                    supplied_comparison_rate=R_DELTA if delta else R,
                    in_supplied_comparison_range=gain >= (R_DELTA if delta else R),
                    assurance='stationary numerical witness; no new interval/global certification',
                    values={key:encode(value) for key,value in values.items()},
                    vectors=vectors, checks={key:encode(value) for key,value in checks.items()},
                    residual_threshold=encode(threshold), numerical_checks_passed=passed)

def compare_precision(a, b):
    if (a['gain'], a['delta']) != (b['gain'], b['delta']):
        raise ValueError('Cannot compare different cases.')
    with mp.workdps(130):
        diffs = {}
        # Require relative agreement even for the tiny in-domain direction errors/gaps.
        for key in a['values']:
            x, y = mp.mpf(a['values'][key]), mp.mpf(b['values'][key])
            diffs[key] = abs(x-y)/max(abs(y), mp.mpf('1e-100'))
        for key in a['vectors']:
            va = mp.matrix([mp.mpf(x) for x in a['vectors'][key]])
            vb = mp.matrix([mp.mpf(x) for x in b['vectors'][key]])
            diffs[key] = mp.norm(va-vb)/max(1,mp.norm(vb))
        passed = all(value < mp.mpf('1e-55') for value in diffs.values())
        if not passed:
            raise ArithmeticError('Precision agreement failed.')
        return dict(relative_differences={key:mp.nstr(value,12) for key,value in diffs.items()},
                    tolerance='1e-55', passed=passed)

def run(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    specs = [(gain,0) for gain in DIRECTION_GAINS]+[(gain,1) for gain in WEIGHTED_GAINS]
    rows, searches, controls = [], [], []
    for gain, delta in specs:
        if gain < R:
            coarse, fine = coarse_search(gain,delta,1), coarse_search(gain,delta,2)
            seed = fine['selected']['fast_time']
            gap = abs(coarse['selected']['height']-fine['selected']['height'])
            if gap > 2e-10:
                raise ArithmeticError('Seed-grid refinement changed the selected height.')
            search = dict(coarse=coarse, refined=fine, height_difference=gap)
        else:
            seed, search = None, dict(method='high-precision asymptotic fast-time seed; no double-precision large-gain sweep')
        a, b = (compute_case(gain, delta, dps, seed) for dps in PRECISIONS)
        comparison = compare_precision(a,b)
        record = dict(gain=gain, delta=delta, at_80_digits=a, at_110_digits=b,
                      precision_comparison=comparison, search=search)
        if delta:
            # Check the unshifted branch independently at the identical gain.
            base_seed = coarse_search(gain,0,1)['selected']['fast_time'] if gain<R else None
            unshifted = compute_case(gain,0,110,base_seed)
            with mp.workdps(120):
                t0 = mp.mpf(unshifted['values']['stationary_time'])
                tw = mp.mpf(b['values']['stationary_time'])
                m = model(gain,1)
                T0 = mp.expm(m['A']*t0)
                s0,u0,_,_,_ = linked_pair(T0,m['u'],m['v'])
                log_deriv_at_unshifted = (u0.T*m['A']*u0)[0]
                separation = tw-t0
                if not (separation>0 and abs(log_deriv_at_unshifted-1)<mp.mpf('1e-80')):
                    raise ArithmeticError('Shifted-time independence check failed.')
                record['shifted_time_check'] = dict(
                    unshifted_time=mp.nstr(t0,100), weighted_time=mp.nstr(tw,100),
                    weighted_minus_unshifted_time=mp.nstr(separation,100),
                    weighted_log_derivative_at_unshifted_time=mp.nstr(log_deriv_at_unshifted,100),
                    weighted_height_at_unshifted_time=mp.nstr(s0,100), passed=True)
            controls.append(unshifted)
        rows.append(record)
        print(f'k={gain}, delta={delta}, h={b["values"]["stationary_height"][:17]}, '
              f'in-domain={b["in_supplied_comparison_range"]}, precision/residual checks PASS', flush=True)
    result = dict(schema_version=1, component='direction_diagnostics',
        method='rate-scaled log-norm stationarity, matrix exponential, linked symmetric-Gram singular pair',
        assurance='numerical observations, stationary evaluations with no new all-time enclosure',
        precisions=list(PRECISIONS), case_count=len(rows), rows=rows,
        unshifted_time_controls=controls,
        supplied_comparison_rates=dict(unshifted=R,shifted_delta_1=R_DELTA),
        all_passed=True)
    (output/'direction_diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')
    for name,delta in [('direction_diagnostics.csv',0),('weighted_envelope_diagnostics.csv',1)]:
        chosen=[row['at_110_digits'] for row in rows if row['delta']==delta]
        fields=['gain','delta','supplied_comparison_rate','in_supplied_comparison_range','assurance']+list(chosen[0]['values'])
        with (output/name).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields)
            writer.writeheader()
            for row in chosen:
                writer.writerow({**{key:row[key] for key in fields[:5]},**row['values']})
    return result
