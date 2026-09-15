"""Integer-only interval replay of finite-gain certificate inequalities.

The enclosing arithmetic uses Python integers on a fixed binary grid. It does
not use mpmath, NumPy, SciPy, floating-point eigensolvers or transcendental
functions. Recorded spectral candidates are proposals and are checked anew.
The matrix-exponential enclosure uses Taylor tails and scaling/squaring.
This is an independent arithmetic implementation, subject to arithmetic-implementation correctness.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt, factorial
from pathlib import Path
import argparse, json, time, hashlib
from .certificate_checks import (verify_record_structure, endpoint_fraction,
                                   strict_margin_record, rational)
BITS=256
SCALE=1<<BITS

def ceildiv(a:int,b:int)->int:
    return -((-a)//b)

class DI:
    """Closed interval [lo/SCALE, hi/SCALE] with outward integer arithmetic."""
    __slots__=('lo','hi')
    def __init__(self,lo:int,hi:int):
        if lo>hi:raise ValueError('reversed interval')
        self.lo,self.hi=lo,hi
    @classmethod
    def scalar(cls,value):
        q=rational(value)
        return cls(q.numerator*SCALE//q.denominator,
                   ceildiv(q.numerator*SCALE,q.denominator))
    @classmethod
    def range(cls,left,right):
        a,b=rational(left),rational(right)
        if a>b:raise ValueError('reversed interval')
        return cls(cls.scalar(a).lo,cls.scalar(b).hi)
    def __neg__(self):return DI(-self.hi,-self.lo)
    def __add__(self,other):
        if not isinstance(other,DI):other=DI.scalar(other)
        return DI(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __sub__(self,other):return self+-as_di(other)
    def __rsub__(self,other):return as_di(other)+-self
    def __mul__(self,other):
        other=as_di(other)
        vals=[self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi]
        return DI(min(vals)//SCALE,ceildiv(max(vals),SCALE))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=as_di(other)
        if other.lo<=0<=other.hi:raise ZeroDivisionError('interval denominator contains zero')
        rec=DI(SCALE*SCALE//other.hi,ceildiv(SCALE*SCALE,other.lo))
        return self*rec
    def __abs__(self):
        if self.lo>=0:return DI(self.lo,self.hi)
        if self.hi<=0:return -self
        return DI(0,max(-self.lo,self.hi))
    def __pow__(self,n):
        if isinstance(n,bool) or not isinstance(n,int) or n<0:raise ValueError('nonnegative integer power')
        ans=DI.scalar(1)
        for _ in range(n):ans=ans*self
        return ans
    def sqrt(self):
        if self.lo<0:raise ValueError('negative square-root interval')
        lo=isqrt(self.lo*SCALE);z=self.hi*SCALE;hi=isqrt(z)
        if hi*hi<z:hi+=1
        return DI(lo,hi)
    def endpoints(self):return F(self.lo,SCALE),F(self.hi,SCALE)

def as_di(x):return x if isinstance(x,DI) else DI.scalar(x)
def mat(values):return [[as_di(x) for x in row] for row in values]
def eye(n):return mat([[int(i==j) for j in range(n)] for i in range(n)])
def trans(A):return list(map(list,zip(*A)))
def add(A,B):return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def scale(A,z):return [[a*z for a in row] for row in A]
def mul(A,B):
    if len(A[0])!=len(B):raise ValueError('matrix dimensions')
    Bt=trans(B)
    return [[sum((x*y for x,y in zip(ar,bc)),DI.scalar(0)) for bc in Bt] for ar in A]

def infinity_upper(A):
    v=max(sum(abs(x).hi for x in row) for row in A)
    return DI(v,v)

def expm(A,degree=34):
    n=len(A)
    if n<1 or any(len(row)!=n for row in A):raise ValueError('nonempty square matrix')
    if isinstance(degree,bool) or not isinstance(degree,int) or degree<6:raise ValueError('degree')
    Z=A;v=infinity_upper(Z);s=0
    while v.hi>SCALE//2:
        Z=scale(Z,F(1,2));v=infinity_upper(Z);s+=1
        if s>2048:raise RuntimeError('excessive scaling')
    term=eye(n);total=eye(n)
    for j in range(1,degree+1):
        term=scale(mul(term,Z),F(1,j));total=add(total,term)
    rho=(v**(degree+1))/factorial(degree+1)/(1-v/(degree+2))
    tail=DI(-rho.hi,rho.hi)
    total=[[z+tail for z in row] for row in total]
    for _ in range(s):total=mul(total,total)
    return total

def positive_definite(A):
    n=len(A)
    if not n or any(len(row)!=n for row in A):raise ValueError('square matrix')
    B=[row[:] for row in A]
    for i in range(n):
        for j in range(i):
            lo=max(B[i][j].lo,B[j][i].lo);hi=min(B[i][j].hi,B[j][i].hi)
            if lo>hi:return False
            B[i][j]=B[j][i]=DI(lo,hi)
    L=eye(n);D=[]
    for i in range(n):
        pivot=B[i][i]-sum((L[i][j]*L[i][j]*D[j] for j in range(i)),DI.scalar(0))
        if pivot.lo<=0:return False
        D.append(pivot)
        for row in range(i+1,n):
            L[row][i]=(B[row][i]-sum((L[row][j]*L[i][j]*D[j] for j in range(i)),DI.scalar(0)))/pivot
    return True

def model(metric):
    if metric=='raw':return mat([[0,1,0],[0,0,0],[1,0,-5]]),mat([[1],[3],[0]]),mat([[-1],[0],[F(1,2)]])
    if metric=='normalized':return mat([[0,5,0],[0,0,0],[F(1,2),0,-5]]),mat([[1],[F(3,5)],[0]]),mat([[-1],[0],[1]])
    raise ValueError('unknown metric')

def unit(values):
    q=[rational(z) for z in values];s=sum(z*z for z in q)
    if s<=0:raise ValueError('zero witness')
    return [[DI.scalar(z)/DI.scalar(s).sqrt()] for z in q]

def decimal_interval(z,digits=18):
    def fmt(n):
        sign='-' if n<0 else '';text=str(abs(n)).zfill(digits+1)
        return sign+(text[:-digits]+'.'+text[-digits:] if digits else text)
    return {'lower':fmt(z.lo*10**digits//SCALE),
            'upper':fmt(ceildiv(z.hi*10**digits,SCALE))}

def replay_case(record,case):
    verify_record_structure(record,case)
    A0,b,c=model(case['metric']);bc=mul(b,trans(c));A=add(A0,scale(bc,F(case['gain'])))
    n=len(A);u_record=record['upper_certificate']
    # All proposed bounds are frozen exact rationals, not accepted as evidence.
    mu=endpoint_fraction(u_record['mu_upper']['binary_endpoints'][1])
    H=scale(add(A,trans(A)),F(1,2))
    if not positive_definite(add(scale(eye(n),mu),scale(H,-1))):raise ValueError('mu proposal rejected')
    best=F(1);grid=[];stop=None
    for row in u_record['samples']:
        j=row['step'];E=expm(scale(A,F(case['step'])*j))
        z=endpoint_fraction(row['norm_upper']['binary_endpoints'][1])
        if z<=0 or not positive_definite(add(scale(eye(n),z*z),scale(mul(trans(E),E),-1))):
            raise ValueError('sample norm proposal rejected')
        best=max(best,z);stop=z
        grid.append({'sample':j,'norm_candidate_verified':True})
    if stop is None or stop>=1:raise ValueError('tail is not strictly contractive')
    factor=expm(mat([[max(F(0),mu)*F(case['step'])]]))[0][0]
    ub=factor*best
    if ub.endpoints()[1]>=F(case['display_upper']):raise ValueError('upper claim not recovered')
    # Recompute every lower-bound cell rather than trust recorded witness values.
    d=-(mul(trans(c),b)[0][0]);Q=[[-x/d for x in row] for row in bc]
    u=unit(case['u_rational']);v=unit(case['v_rational']);K=DI.scalar(case['cutoff']);tau=F(case['tau'])
    cells=[];lower=None
    for j in range(case['cells']):
        theta=DI.range(F(j,case['cells']),F(j+1,case['cells']))
        M=scale(add(scale(Q,-1),scale(A0,theta/(d*K))),tau)
        z=mul(mul(trans(u),expm(M)),v)[0][0]
        lower=z.lo if lower is None else min(lower,z.lo)
        cells.append({'cell':j,'x_left':str(F(j,case['cells'])),'x_right':str(F(j+1,case['cells'])),
                      'witness_enclosure':decimal_interval(z)})
    lb=DI(lower,lower)
    if lb.endpoints()[0]<=F(case['display_lower']):raise ValueError('lower claim not recovered')
    return {'name':case['name'],'arithmetic':'integer-only fixed-grid dyadic intervals',
            'binary_bits':BITS,'taylor_degree':34,'upper':decimal_interval(ub),'lower':decimal_interval(lb),
            'positive_time_samples':len(grid),'time_zero_included':True,
            'tail_norm_verified_below_one':True,'grid_checks':grid,'cells':cells,
            'strict_margin':strict_margin_record(case['gain'],case['eta']),
            'all_public_inequalities_verified':True,
            'formal_software_verification':False}

