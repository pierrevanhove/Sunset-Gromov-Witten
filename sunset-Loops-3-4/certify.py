"""Rigorous real-ball continuation for four-loop endpoint connections.
Every truncation has an explicit majorant; no floating-point input data.
"""
from flint import arb,arb_mat,ctx
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import json,time,argparse
Q=[[0,0,0,0,1],[5,28,63,70,35],[285,1088,1580,1036,259],[900,2700,2925,1350,225]]
S=[[int(k==r==0) for r in range(5)] for k in range(5)]
for k in range(1,5):
    for r in range(1,k+1):S[k][r]=S[k-1][r-1]+r*S[k-1][r]
ODE=[[sum(Q[j][k]*S[k][r] for k in range(r,5) for j in range(4) if d==j+r) for d in range(r+4)] for r in range(5)]
def A(v):return arb(v.numerator)/v.denominator if isinstance(v,F) else arb(v)
def zero(n):return [arb(0) for _ in range(n)]
def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,c):return [x*c for x in a]
def deriv(p):return [(k+1)*p[k+1] for k in range(len(p)-1)]+[arb(0)]
def op(p,k,r):
    for _ in range(r):p=add(scale(p,k),deriv(p))
    return p

def pv(p,x):
    v=arb(0)
    for a in reversed(p):v=v*x+a
    return v

def mul(a,b):return [sum((a[j]*b[k-j] for j in range(k+1)),arb(0)) for k in range(5)]
def qshift(q,d):return [sum(q[a]*comb(a,k)*d**(a-k) for a in range(k,5)) for k in range(5)]
def norm(p):return sum((x.abs_upper() for x in p),arb(0)).upper()
def maxabs(p):return max(x.abs_upper() for x in p)
def widen(x,e):return x+arb(0,e.abs_upper())

def germ_coeff(N):
    us=[[arb(1)]+zero(4)]
    for d in range(1,N+1):
        rhs=zero(5)
        for l in range(1,min(3,d)+1):rhs=add(rhs,mul(list(map(arb,qshift(Q[l],d-l))),us[d-l]))
        inv=[arb((-1)**r*comb(r+3,3))/arb(d)**(4+r) for r in range(5)]
        us.append(scale(mul(inv,rhs),-1))
    # Base cases of the coefficient majorant. Strict bound 2*100^d allows rounding.
    for d in range(10):assert norm(us[d]) < 2*arb(100)**d
    assert sum(F(sum(Q[l][a]*10**a for a in range(5)),9**4*100**l) for l in range(1,4))<1
    return us

def initial(N):
    y=arb(1)/512;t=y.log();us=germ_coeff(N);W=[];C=[]
    q=arb(100)/512
    # Uniform tail bound for theta jets of order <=3 and primitive.
    e=1024*q**(N+1)*(N+2)**3*(1+4*q+q*q)/(1-q)**4
    for j in range(5):
        jets=zero(4);c=arb(0)
        for d,row in enumerate(us):
            p=[row[j-r]/factorial(r) for r in range(j+1)]
            for r in range(4):jets[r]+=y**d*pv(op(p,d,r),t)
            if d:
                prim=zero(len(p));dp=p
                for r in range(len(p)):
                    prim=add(prim,scale(dp,arb((-1)**r)/arb(d)**(r+1)));dp=deriv(dp)
                c+=y**d*pv(prim,t)
            else:c+=t**(j+1)/factorial(j+1)
        jets=[widen(v,e) for v in jets];c=widen(c,e)
        W.append([jets[0],jets[1]/y,(jets[2]-jets[1])/y**2,(jets[3]-3*jets[2]+2*jets[1])/y**3]);C.append(c)
    return y,W,C,e

def shifted(poly,y):return [sum((arb(poly[j])*comb(j,k)*y**(j-k) for j in range(k,len(poly))),arb(0)) for k in range(len(poly))]
def step(y,W,C,h,N):
    aa=[shifted(p,y) for p in ODE];r=arb(1)/4;q=arb(1)/4
    # Supremum of augmented logarithmic-jet first-order system on |z|<=r.
    ymax=y*(1+r);den=(1+y*(1-r))*(1+9*y*(1-r))*(1+25*y*(1-r))
    numer=arb(1)+sum((sum(Q[l][:4])*ymax**l for l in range(1,4)),arb(0))
    M=max(arb(1),(numer/den).upper())/(1-r)
    nw=[];nc=[];largest=arb(0)
    for j,jet in enumerate(W):
        theta=[jet[0],y*jet[1],y*y*jet[2]+y*jet[1],y**3*jet[3]+3*y*y*jet[2]+y*jet[1]]
        B=maxabs(theta+[C[j],arb(1)])*(M*r).exp()
        e=B*q**(N+1)*(1+4*q+q*q)/(1-q)**4
        largest=max(largest,e.upper())
        f=[jet[k]/factorial(k) for k in range(4)]
        for n in range(N-3):
            rhs=arb(int(j==4 and n==0))
            for rr in range(5):
                for k,a in enumerate(aa[rr]):
                    if k>n or rr==4 and k==0:continue
                    ix=n-k+rr
                    factor=1
                    for v in range(n-k+1,ix+1):factor*=v
                    rhs-=a*f[ix]*factor
            f.append(rhs/(aa[4][0]*(n+1)*(n+2)*(n+3)*(n+4)))
        vals=[];df=f
        for rr in range(4):
            vals.append(widen(pv(df,h),e*(4*(N+1)/(r*y))**rr));df=deriv(df)
        cf=[C[j]];prev=arb(0)
        for n in range(N):
            qq=(f[n]-prev)/y;cf.append(qq/(n+1));prev=qq
        nw.append(vals);nc.append(widen(pv(cf,h),e))
    return nw,nc,largest

def frob(v,source,N):
    ps=[list(map(arb,v[:2])),list(map(arb,v[2:]))]
    for k in range(2,N+1):
        p=ps[k-1];rhs=add(scale(op(p,k-1,4),259),scale(op(p,k-1,2),26))
        for a,c in enumerate(Q[1]):rhs=add(rhs,scale(op(ps[k-2],k-2,a),c))
        if k>=3:rhs=add(rhs,op(ps[k-3],k-2,4))
        if k==2:rhs[0]-=source
        rhs=add(rhs,scale(deriv(rhs),-arb(2)/k-arb(2)/(k-1)))
        ps.append(scale(rhs,-arb(1)/(225*k*k*(k-1)**2)))
    for k in range(4):assert norm(ps[k]) < 2*arb(8)**k
    return ps

def eval_frob(ps,x):
    N=len(ps)-1;t=x.log();q=8*x;jets=zero(4);tail=arb(0)
    for k,p in enumerate(ps):
        for r in range(4):jets[r]+=x**(k+1)*pv(op(p,k+1,r),t)
        tail+=x**(k+1)*pv(add(p,scale(deriv(p),-arb(1)/(k+1))),t)/(k+1)
    logbound=max(arb(1),t.abs_upper())
    base=2*x*logbound*q**(N+1)*(1+4*q+q*q)/(1-q)**4
    jets=[widen(v,base*(N+3)**r) for r,v in enumerate(jets)]
    tail=widen(tail,4*x*logbound*q**(N+1)/(1-q))
    return jets,tail

def box(v):
    den=10**24
    lo=int((v.lower()*den).lower().floor().unique_fmpz())
    hi=int((v.upper()*den).upper().ceil().unique_fmpz())
    assert v >= arb(lo)/den and v <= arb(hi)/den
    return {'ball':v.str(45),'outer_lower_numerator':str(lo),'outer_upper_numerator':str(hi),'outer_denominator':str(den),'radius':v.rad().str(8)}
def run(bits=1400,N=400,Ng=450,Nf=250,match=32):
    assert match>=32 and N>=4 and Ng>=10 and Nf>=4
    ctx.prec=bits;start=time.time();y,W,C,initial_tail=initial(Ng);steps=0;maxrem=arb(0)
    while y<match:
        h=min(y/16,arb(match)-y);W,C,e=step(y,W,C,h,N);y+=h;steps+=1;maxrem=max(maxrem,e)
        if steps%40==0:print('step',steps,'y',y.str(8),'radius',maxabs([a.rad() for w in W for a in w]).str(5),flush=True)
    x=1/y;bases=[frob([int(a==j) for a in range(4)],0,Nf) for j in range(4)];part=frob([0]*4,1,Nf)
    be=[eval_frob(b,x) for b in bases];pe=eval_frob(part,x)
    matrix=arb_mat([[be[j][0][r] for j in range(4)] for r in range(4)])
    connections=[];limits=[]
    for j,w in enumerate(W):
        target=[w[0],-y*w[1],y*w[1]+y*y*w[2],-y*w[1]-3*y*y*w[2]-y**3*w[3]]
        if j==4:target=[a-b for a,b in zip(target,pe[0])]
        sol=matrix.solve(arb_mat([[v] for v in target]));weights=[sol[k,0] for k in range(4)]
        connections.append(weights);limits.append(C[j]+sum((weights[k]*be[k][1] for k in range(4)),arb(0))+(pe[1] if j==4 else 0))
    R=limits[0];a0,b0=connections[0][:2]
    assert b0<0
    ds=[(v[0]*b0-a0*v[1])/b0**2 for v in connections]
    Es=[sum(((-R)**(k-1-h)/factorial(k-1-h)*ds[h] for h in range(1,k)),arb(0)) for k in range(2,6)]
    for e,sign in zip(Es,[1,-1,1,1]):assert sign*e>0
    out={'bits':bits,'taylor_order':N,'germ_order':Ng,'frobenius_order':Nf,'matching_y':match,'steps':steps,'initial_tail_bound':box(initial_tail),'largest_step_tail_W_bound':box(maxrem),'R_infinity':box(R),'Q_infinity':box(R.exp()),'connections':[[box(v) for v in row] for row in connections],'C_limits':[box(v) for v in limits],'E':[box(v) for v in Es],'seconds':time.time()-start,'certificate_passed':True}
    Path('certificate_'+str(match)+'.json').write_text(json.dumps(out,indent=2)+'\n')
    print('R',R.str(40),'E',[v.str(30) for v in Es],'seconds',out['seconds'],flush=True)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bits',type=int,default=1400);p.add_argument('--order',type=int,default=400);p.add_argument('--germ',type=int,default=450);p.add_argument('--frob',type=int,default=250);p.add_argument('--match',type=int,default=32);a=p.parse_args();run(a.bits,a.order,a.germ,a.frob,a.match)

