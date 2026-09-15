"""Small-matrix interval certificates for continuous-time transient peaks.

Certificate operations use mpmath.iv directed interval arithmetic. Floating
point proposes spectral bounds and witness vectors; interval LDL or interval
normalization validates them. Matrix exponentials use an explicit Taylor tail
and interval scaling/squaring. Assurance is conditional on the arithmetic
implementation and correct execution, not a formal software audit."""
from __future__ import annotations
import math
from fractions import Fraction
import mpmath as mp
import numpy as np
from scipy.linalg import expm,svd
iv=mp.iv

def precision(dps:int=50)->None:
    if isinstance(dps,bool) or not isinstance(dps,int) or dps<30:
        raise ValueError('Use an integer of at least 30 decimal interval digits.')
    iv.dps=dps
    mp.mp.dps=dps+35  # Reporting/proposals only. Never used for enclosure arithmetic.
precision()

def exact(value):
    """Enclose an exactly specified finite rational. Reject binary floats."""
    if isinstance(value,bool) or not isinstance(value,(int,str,Fraction)):
        raise TypeError('Use an integer, rational string, or Fraction, not a float.')
    rational=Fraction(value)  # Rejects infinities, NaNs and interval literals.
    return iv.mpf(rational.numerator)/rational.denominator

def bounds(x)->tuple[mp.mpf,mp.mpf]:
    # Higher mp precision represents finite iv endpoint tuples exactly.
    return mp.mpf(x._mpi_[0]),mp.mpf(x._mpi_[1])

def endpoint_fraction(endpoint)->Fraction:
    """Decode a finite binary endpoint without floating-point conversion."""
    sign,mantissa,exponent,bitcount=map(int,endpoint)
    if sign not in (0,1) or mantissa<0 or bitcount<0:
        raise ValueError('A finite binary endpoint is required.')
    val=Fraction((-1 if sign else 1)*mantissa)
    return val*(2**exponent) if exponent>=0 else val/Fraction(2**(-exponent))

def _decimal_grid(integer:int,digits:int)->str:
    sign='-' if integer<0 else ''
    text=str(abs(integer)).zfill(digits+1)
    return sign+(text[:-digits]+'.'+text[-digits:] if digits else text)

def display_bounds(x,digits:int=18)->dict:
    """Outward decimal serialization using exact integer floor and ceiling."""
    if isinstance(digits,bool) or not isinstance(digits,int) or digits<0:
        raise ValueError('digits must be a nonnegative integer')
    aa,bb=(endpoint_fraction(t) for t in x._mpi_)
    scale=10**digits
    lo=(aa.numerator*scale)//aa.denominator
    hi=-((-bb.numerator*scale)//bb.denominator)
    return {'lower':_decimal_grid(lo,digits),'upper':_decimal_grid(hi,digits),
            'binary_endpoints':[[int(y) for y in z] for z in x._mpi_]}

def norminf_upper(A):
    out=iv.mpf(0)
    for i in range(A.rows):
        row=sum((abs(A[i,j]) for j in range(A.cols)),iv.mpf(0)).b
        if row>out:out=row
    return out

def interval_expm(A,degree:int=26):
    if A.rows!=A.cols:raise ValueError('square matrix required')
    if isinstance(degree,bool) or not isinstance(degree,int) or degree<6:
        raise ValueError('Taylor degree must be an integer of at least six')
    if A.rows<1:raise ValueError('nonempty matrix required')
    n=A.rows;Z=A.copy();v=norminf_upper(Z);scale=0
    while v>iv.mpf('0.5'):
        Z=Z/2;v=(v/2).b;scale+=1
        if scale>2048:raise RuntimeError('Excessive scaling, certificate refused')
    result=iv.eye(n);term=iv.eye(n)
    for j in range(1,degree+1):
        term=term*Z/j;result=result+term
    rho=(v**(degree+1)/math.factorial(degree+1)/(1-v/(degree+2))).b
    # ||Taylor tail||_infinity <= rho, hence every entry lies in [-rho,rho].
    slack=iv.mpf([-rho,rho])
    for i in range(n):
        for j in range(n):result[i,j]+=slack
    for _ in range(scale):result=result*result
    return result

def mid_float(A)->np.ndarray:
    return np.array([[float(sum(bounds(A[i,j]))/2) for j in range(A.cols)] for i in range(A.rows)])

def positive_definite(B)->bool:
    """Sufficient interval LDL check for every symmetric matrix enclosed by B."""
    n=B.rows
    if n<1 or n!=B.cols:raise ValueError('nonempty square matrix required')
    # Only certify a symmetric target. Intersect paired off-diagonal ranges.
    B=B.copy()
    for i in range(n):
        for j in range(i):
            lo=max(B[i,j].a,B[j,i].a);hi=min(B[i,j].b,B[j,i].b)
            if lo>hi:return False
            B[i,j]=B[j,i]=iv.mpf([lo,hi])
    L=iv.eye(n);D=[]
    for i in range(n):
        pivot=B[i,i]-sum((L[i,j]*L[i,j]*D[j] for j in range(i)),iv.mpf(0))
        if not pivot.a>0:return False
        D.append(pivot)
        for r in range(i+1,n):
            L[r,i]=(B[r,i]-sum((L[r,j]*L[i,j]*D[j] for j in range(i)),iv.mpf(0)))/pivot
    return True

def spectral_upper(E):
    """Validated norm upper bound. Midpoint SVD only proposes the candidate."""
    proposal=float(np.linalg.norm(mid_float(E),2));pad=max(1e-13,abs(proposal)*1e-12)
    EE=E.T*E
    for _ in range(80):
        z=exact(format(proposal+pad,'.17g'))
        if positive_definite(z*z*iv.eye(E.cols)-EE):return z
        pad*=2
    raise RuntimeError('Could not verify a finite spectral upper bound')

def mu_upper(A):
    H=(A+A.T)/2;proposal=float(np.linalg.eigvalsh(mid_float(H))[-1]);pad=max(1e-13,abs(proposal)*1e-12)
    for _ in range(80):
        z=exact(format(proposal+pad,'.17g'))
        if positive_definite(z*iv.eye(A.rows)-H):return z
        pad*=2
    raise RuntimeError('Could not certify numerical-abscissa upper bound')

def exact_pll(metric:str='raw'):
    A0=iv.matrix([[0,1,0],[0,0,0],[1,0,-5]])
    b=iv.matrix([1,3,0]);c=iv.matrix([-1,0,exact('0.5')])
    if metric=='normalized':
        S=iv.diag([exact(1),exact(1)/5,exact(1)/2]);Si=iv.diag([1,5,2])
        A0=S*A0*Si;b=S*b;c=Si.T*c
    elif metric!='raw':raise ValueError('unknown metric')
    return A0,b,c

def fixed_gain_upper(k:str,step:str,metric:str='raw',max_steps:int=5000)->dict:
    """Verified continuous-time upper bound, no low-rank truncation.

    Each sample is independently interval-enclosed. A positive norm<1 stopping
    sample establishes the infinite-tail reduction. No sampled time is skipped.
    """
    if isinstance(max_steps,bool) or not isinstance(max_steps,int) or max_steps<1:
        raise ValueError('max_steps must be a positive integer')
    A0,b,c=exact_pll(metric);K=exact(k);h=exact(step)
    if not (K.a>0 and h.a>0):raise ValueError('positive gain and step required')
    A=A0+K*b*c.T;mu=mu_upper(A)
    best=exact(1);rows=[]
    if mu.b<=0:return {'gain':k,'metric':metric,'upper':display_bounds(best),'steps':0,'contraction':True}
    for j in range(1,max_steps+1):
        E=interval_expm((j*h)*A);ub=spectral_upper(E)
        if ub.b>best:best=ub.b
        rows.append({'step':j,'norm_upper':display_bounds(ub,14)})
        if ub.b<1:
            total=(best*iv.exp(mu*h)).b
            return {'gain':k,'metric':metric,'step':step,'steps':j,'stop_time':display_bounds(j*h),
                    'mu_upper':display_bounds(mu),'grid_upper':display_bounds(best),'upper':display_bounds(total),
                    'stop_norm_upper':display_bounds(ub),'samples':rows,'interval_dps':iv.dps,
                    'taylor_degree':26,'includes_time_zero':True,
                    'status':'directed-interval certificate; Taylor tail plus interval LDL plus semigroup tail'}
    raise RuntimeError('Stopping condition not verified; no certificate returned')

def normalized_rational_vector(values:list[str]):
    if not isinstance(values,(list,tuple)) or not values:
        raise ValueError('a nonempty vector is required')
    v=iv.matrix([exact(x) for x in values]);s=(v.T*v)[0]
    if not s.a>0:raise ValueError('nonzero vector required')
    return v/iv.sqrt(s)

def high_gain_witness(K:str,tau:str,metric:str='raw',cells:int=128,uv:tuple[list[str],list[str]]|None=None)->dict:
    """Certify a lower bound for every k>=K by x=K/k in [0,1].

    t= tau/(d*k) is an admissible witness time for every finite k.
    Cells are rational closed intervals; the limiting point x=0 is included.
    No stability or monotonicity of the actual peak is assumed.
    """
    if isinstance(cells,bool) or not isinstance(cells,int) or cells<1:
        raise ValueError('positive integer cell count required')
    A0,b,c=exact_pll(metric);Kiv=exact(K);tauiv=exact(tau);d=-(c.T*b)[0]
    if not (Kiv.a>0 and tauiv.a>=0 and d.a>0):raise ValueError('invalid witness parameters')
    Q=-(b*c.T)/d
    if uv is None:
        T=expm(mid_float((-Q+A0/(d*Kiv))*tauiv));U,s,Vh=svd(T)
        uv=([format(z,'.16g') for z in U[:,0]],[format(z,'.16g') for z in Vh[0,:]])
    if len(uv)!=2 or any(len(w)!=A0.rows for w in uv):
        raise ValueError('witness dimensions must equal the state dimension')
    u=normalized_rational_vector(uv[0]);v=normalized_rational_vector(uv[1]);rows=[];lo=None
    for j in range(cells):
        xa=exact(j)/cells;xb=exact(j+1)/cells
        x=iv.mpf([xa.a,xb.b]);H=(-Q+x*A0/(d*Kiv))*tauiv
        T=interval_expm(H);value=(u.T*T*v)[0]
        if lo is None or value.a<lo:lo=value.a
        rows.append({'cell':j,'x_left':f'{j}/{cells}','x_right':f'{j+1}/{cells}',
                     'scalar_witness':display_bounds(value,16)})
    return {'cutoff':K,'tau':tau,'metric':metric,'cells':cells,'u_rational':uv[0],'v_rational':uv[1],
            'lower':display_bounds(lo),'interval_dps':iv.dps,'cell_records':rows,
            'taylor_degree':26,'normalization':'exact rational vectors divided by their exact Euclidean norms',
            'status':'directed-interval lower certificate for the entire half-line k>=K'}
