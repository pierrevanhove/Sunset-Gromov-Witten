"""Exact normalization checks and rational checks on exported enclosures."""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import json
from certify import Q

def conv(a,b,N):
    out=[F(0)]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:N+1-i]):out[i+j]+=x*y
    return out

def exact_recurrence(N):
    out=[[F(1)]+[F(0)]*4]
    for d in range(1,N+1):
        rhs=[F(0)]*5
        for l in range(1,min(d,3)+1):
            p=[F(sum(Q[l][a]*comb(a,k)*(d-l)**(a-k) for a in range(k,5))) for k in range(5)]
            v=conv(p,out[d-l],4);rhs=[a+b for a,b in zip(rhs,v)]
        inv=[F((-1)**k*comb(k+3,3),d**(k+4)) for k in range(5)]
        out.append([-v for v in conv(rhs,inv,4)])
    return out

def squarefree(N):
    f=[F(1,factorial(d)**2) for d in range(N+1)];g=[];h=F(0)
    for d in range(N+1):
        if d:h+=F(1,d)
        g.append(-2*h*f[d])
    ps=[]
    for k in range(5):
        a=[F(1)]+[F(0)]*N
        for i in range(5):a=conv(a,g if i<k else f,N)
        ps.append(a)
    out=[];num=[F(1)]+[F(0)]*4
    for d in range(N+1):
        if d:
            for _ in range(2):num=[num[k]+(num[k-1]/d if k else 0) for k in range(5)]
        out.append([(-1)**d*factorial(d)**2*sum((num[a]*ps[k-a][d]/factorial(k-a) for a in range(k+1)),F(0)) for k in range(5)])
    return out

def interval(b):return F(int(b['outer_lower_numerator']),int(b['outer_denominator'])),F(int(b['outer_upper_numerator']),int(b['outer_denominator']))

if __name__=='__main__':
    a=exact_recurrence(9);b=squarefree(9);assert a==b
    base=[sum(map(abs,row)) <=2*100**d for d,row in enumerate(a)];assert all(base)
    germ_factor=sum(F(sum(Q[l][k]*10**k for k in range(5)),9**4*100**l) for l in range(1,4));assert germ_factor<1
    H1=(F(259)+F(26,16))*F(16,9)*F(13,6)/225
    H2=F(201,225)*F(13,6);H3=F(1,225)*F(13,6)
    assert H1<5 and H2<2 and H3<F(1,100)
    certs=[json.loads(Path('certificate_'+str(m)+'.json').read_text()) for m in [32,48]]
    Es=[]
    for c in certs:
        assert c['certificate_passed']
        bounds=[interval(v) for v in c['E']]
        assert bounds[0][0]>0 and bounds[1][1]<0 and bounds[2][0]>0 and bounds[3][0]>0
        Es.append(bounds)
    for left,right in zip(Es[0],Es[1]):assert max(left[0],right[0])<=min(left[1],right[1])
    out={'exact_squarefree_normalization_coefficients':50,'germ_base_cases':10,'germ_induction_factor':str(germ_factor),'frobenius_H1_bound':str(H1),'frobenius_H2_bound':str(H2),'frobenius_H3_bound':str(H3),'two_certificates_pass':True,'exported_rational_E_intervals':[[str(a),str(b)] for a,b in Es[0]],'status':'All assertions passed. The analytic tail proofs in the report are part of the certificate.'}
    Path('validation.json').write_text(json.dumps(out,indent=2)+'\n')
    Path('exact_germ_check.json').write_text(json.dumps([[str(v) for v in row] for row in a],indent=2)+'\n')
    print(json.dumps(out,indent=2))

