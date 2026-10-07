"""
Uniform sunset family:  F_L = (1,...,1) hypersurface in (P^1)^{L+1},  X_L = Tot(K_{F_L}).

  L=2 : F_2 = sextic del Pezzo = Bl_3 P^2   (toric)      X_2 local CY 3-fold
  L=3 : F_3 = Mori-Mukai 4-1 threefold      (NOT toric)  X_3 local CY 4-fold
  L=4 : F_4 = Fano fourfold in (P^1)^5      (NOT toric)  X_4 local CY 5-fold

Nothing is hardcoded per L: the Picard-Fuchs operator P_y = sum_j y^j P_j(theta) is
recovered by exact linear algebra from the sunset numbers, so the script extends to any L.

gamma_L = H^{L-1}/(L-1)! = e_{L-1}(p_1,...,p_{L+1})   (integral for every L)
F_gamma = -((L+1)!/(L-1)!) K,  K = z^{-2} coefficient of e^{-HT/z} I^0_X
Identity:  f_{L-1} = (-1)^L theta_Q F_gamma.
"""
from fractions import Fraction as Fr
from math import factorial

N = 16   # order of the Q-series
M = 70   # how many sunset numbers to generate (for finding the operator)

# ================= truncated power series over Q (length N) =================
def mul(a, b):
    r = [Fr(0)]*N
    for i, ai in enumerate(a):
        if ai:
            for j in range(0, N-i):
                if b[j]:
                    r[i+j] += ai*b[j]
    return r

def inv(a):
    assert a[0] != 0
    r = [Fr(0)]*N; r[0] = 1/a[0]
    for n in range(1, N):
        s = sum(a[k]*r[n-k] for k in range(1, n+1))
        r[n] = -s/a[0]
    return r

def compose(f, g):
    assert g[0] == 0
    r = [Fr(0)]*N; p = [Fr(0)]*N; p[0] = Fr(1)
    for n in range(N):
        if f[n]:
            for i in range(N):
                r[i] += f[n]*p[i]
        p = mul(p, g)
    return r

def revert(g):
    h = [Fr(0)]*N; h[1] = 1/g[1]
    for n in range(2, N):
        h[n] -= compose(g, h)[n]/g[1]
    return h

def expser(a):
    assert a[0] == 0
    r = [Fr(0)]*N; r[0] = Fr(1); t = [Fr(0)]*N; t[0] = Fr(1)
    for k in range(1, N):
        t = mul(t, a)
        f = Fr(1, factorial(k))
        for i in range(N):
            r[i] += f*t[i]
    return r

def theta(a):
    return [Fr(n)*a[n] for n in range(N)]

HN = [Fr(0)]*(M+5)
for n in range(1, M+5):
    HN[n] = HN[n-1] + Fr(1, n)

# ================= sunset numbers via generating functions =================
def sunset_data(L, upto):
    """A_D = sum_{|d|=D} multinomial(D;d)^2  over L+1 parts ;
       tilde_beta_D = D!(D-1)! sum_{|d|=D} (sum_i h_{d_i}) / prod (d_i!)^2 ."""
    K = upto+1
    def smul(a, b):
        r = [Fr(0)]*K
        for i, ai in enumerate(a):
            if ai:
                for j in range(0, K-i):
                    if b[j]:
                        r[i+j] += ai*b[j]
        return r
    g = [Fr(1, factorial(d)**2) for d in range(K)]
    h = [HN[d]*Fr(1, factorial(d)**2) for d in range(K)]
    gp = [Fr(1)]+[Fr(0)]*(K-1)
    for _ in range(L+1):
        gp = smul(gp, g)
    gL = [Fr(1)]+[Fr(0)]*(K-1)
    for _ in range(L):
        gL = smul(gL, g)
    hb = smul(h, gL)
    A = [int(Fr(factorial(D))**2*gp[D]) for D in range(K)]
    TB = [Fr(0)]*K
    for D in range(1, K):
        TB[D] = Fr(factorial(D)*factorial(D-1))*(L+1)*hb[D]
    return A, TB

# ================= exact linear solve =================
def solve(rows, rhs):
    m, n = len(rows), len(rows[0])
    Aug = [list(r)+[b] for r, b in zip(rows, rhs)]
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if Aug[i][c] != 0), None)
        if p is None:
            continue
        Aug[r], Aug[p] = Aug[p], Aug[r]
        pv = Aug[r][c]
        Aug[r] = [x/pv for x in Aug[r]]
        for i in range(m):
            if i != r and Aug[i][c] != 0:
                f = Aug[i][c]
                Aug[i] = [x - f*y for x, y in zip(Aug[i], Aug[r])]
        piv.append(c); r += 1
        if r == m:
            break
    for i in range(r, m):
        if all(x == 0 for x in Aug[i][:n]) and Aug[i][n] != 0:
            return None
    if len(piv) < n:
        return None
    sol = [Fr(0)]*n
    for i, c in enumerate(piv):
        sol[c] = Aug[i][n]
    return sol

def find_operator(L, a):
    """P_y = theta^L + sum_{j>=1} y^j P_j(theta),  deg P_j <= L.

    The y-degree J is not free: the leading symbol sum_j (lead P_j) y^j vanishes
    exactly at the finite singular points s = (L+1-2i)^2, i = 0..floor((L+1)/2).
    Dropping s=0 (which sits at y=infinity and contributes no factor) leaves
        J = floor(L/2) + 1
    factors, so J = 2,2,3,3 for L = 2,3,4,5.  We search one past that as a guard:
    if J is too small the linear system is inconsistent, if too large it is
    underdetermined (left-multiply the true operator by y), so solve() returns
    None in both cases and only the minimal J survives."""
    Jexp = L//2 + 1
    for J in range(1, Jexp+2):
        nun = J*(L+1)
        rows, rhs = [], []
        for n in range(J, min(len(a), J+2*nun+15)):
            rows.append([Fr((n-j)**k * a[n-j]) for j in range(1, J+1) for k in range(L+1)])
            rhs.append(Fr(-(n**L)*a[n]))
        if len(rows) < nun:
            continue
        sol = solve(rows, rhs)
        if sol is None:
            continue
        P = [[Fr(0)]*(L+1) for _ in range(J+1)]
        P[0][L] = Fr(1)
        for j in range(1, J+1):
            for k in range(L+1):
                P[j][k] = sol[(j-1)*(L+1)+k]
        assert J == Jexp, f"y-degree {J} != predicted floor(L/2)+1 = {Jexp}"
        for n in range(J, len(a)):               # verify on ALL available n
            tot = sum(sum(P[j][k]*(n-j)**k for k in range(L+1))*a[n-j]
                      for j in range(J+1) if n-j >= 0)
            assert tot == 0, (L, J, n)
        return J, P
    raise RuntimeError("no operator found")

def polyfmt(P):
    t = []
    for k in range(len(P)-1, -1, -1):
        if P[k]:
            t.append(f"{P[k]}" + ("" if k == 0 else ("*th" if k == 1 else f"*th^{k}")))
    return (" + ".join(t) or "0").replace("+ -", "- ")

# ================= jets in rho (Taylor coefficients, order 2) =================
def jm(u, v):
    return (u[0]*v[0], u[0]*v[1]+u[1]*v[0], u[0]*v[2]+u[1]*v[1]+u[2]*v[0])
def jpow(u, k):
    r = (Fr(1), Fr(0), Fr(0))
    for _ in range(k):
        r = jm(r, u)
    return r
def jinv(u):
    a, b, c = u
    return (1/a, -b/a**2, b*b/a**3 - c/a**2)
def jpoly(P, u):
    r = (Fr(0), Fr(0), Fr(0)); p = (Fr(1), Fr(0), Fr(0))
    for k in range(len(P)):
        if P[k]:
            r = tuple(x + P[k]*y for x, y in zip(r, p))
        p = jm(p, u)
    return r

def frobenius(L, J, P):
    C = [(Fr(1), Fr(0), Fr(0))]
    for n in range(1, N):
        m = (Fr(n), Fr(1), Fr(0))
        num = (Fr(0), Fr(0), Fr(0))
        for j in range(1, J+1):
            if n-j >= 0:
                mj = (Fr(n-j), Fr(1), Fr(0))
                num = tuple(x - y for x, y in zip(num, jm(jpoly(P[j], mj), C[n-j])))
        C.append(jm(num, jinv(jpow(m, L))))
    return [c[0] for c in C], [c[1] for c in C]

# ================= A-model side =================
def amodel(L):
    A, TB = sunset_data(L, N)
    S = [Fr(0)]*N; U = [Fr(0)]*N
    for D in range(1, N):
        sg = (-1)**D
        S[D] = Fr(sg*A[D], D)
        U[D] = sg*(Fr(A[D], D)*(HN[D]+HN[D-1]) - Fr(2, L+1)*TB[D])
    SS = mul(S, S)
    K = [U[i] - SS[i]/2 for i in range(N)]
    Q = mul([Fr(0), Fr(1)]+[Fr(0)]*(N-2), expser(S))
    return S, K, Q, revert(Q)

# ================= report =================
NAME = {2: "sextic del Pezzo Bl_3 P^2 (toric)",
        3: "Mori-Mukai 4-1 threefold  (not toric)",
        4: "Fano fourfold in (P^1)^5  (not toric)"}

def fano_name(L):
    """F_L is the (1,...,1) hypersurface in (P^1)^{L+1} for every L."""
    return NAME.get(L, f"(1,...,1) hypersurface in (P^1)^{L+1}")

def report(L):
    print("="*76)
    print(f"  L = {L}   F_L = {fano_name(L)}   ->  X_L local CY {L+1}-fold")
    print("="*76)
    Abig, _ = sunset_data(L, M)
    a = [(-1)**n*Abig[n] for n in range(M+1)]
    print(f"  sunset numbers A_n : {Abig[:7]}")
    J, P = find_operator(L, a)
    print(f"  Picard-Fuchs operator: J = {J}   ({J+1}-term recurrence)")
    for j in range(J+1):
        print(f"     P_{j}(th) = {polyfmt(P[j])}")
    lead = [P[j][L] for j in range(J+1)]
    import numpy as np
    rts = np.roots([float(x) for x in lead][::-1])
    print(f"  finite singular points s = -1/y : "
          f"{sorted(round(float(-1/r),4) for r in rts)}    expect {[(L+1-2*i)**2 for i in range(J)]}")

    RHS = (-1)**(L+1)*factorial(L+1)
    print(f"  inhomogeneity (-1)^(L+1)(L+1)! = {RHS}")

    S, K, Q, yQ = amodel(L)
    KQ = compose(K, yQ)
    cn = Fr(factorial(L+1), factorial(L-1))
    Fg = [-cn*x for x in KQ]
    print(f"\n  K(Q)     = {[str(x) for x in KQ[1:6]]}")
    print(f"  F_gamma  = -{cn}*K = {[str(x) for x in Fg[1:6]]}")

    A0, B = frobenius(L, J, P)
    r1 = mul(B, inv(A0))
    lc = Fr(RHS, factorial(L))
    f = compose([lc*L*(r1[i]-S[i]) for i in range(N)], yQ)
    print(f"  f_{L-1}(Q)  = {[str(x) for x in f[1:6]]}")
    sgn = (-1)**L
    ok = all(f[i] == sgn*theta(Fg)[i] for i in range(1, N))
    print(f"\n  >>> f_{L-1} == {'+' if sgn>0 else '-'}theta_Q F_gamma ?   {ok}    (exact, to Q^{N-1})")

    print(f"\n   m   N_m(gamma_{L})")
    for m in range(1, 9):
        print(f"  {m:2d}   {str(Fg[m]):>18s}")

    # multi-cover tests.  L=2: N^loc = N(gamma)/m by the divisor axiom.
    base = [Fr(0)] + [Fg[d]/d if L == 2 else Fg[d] for d in range(1, N)]
    print("\n  multi-cover test   base_d = sum_{k|d} n_{d/k}/k^w :")
    for w in (1, 2, 3, 4):
        n = {}; good = True
        for d in range(1, N):
            v = base[d]
            for k in range(2, d+1):
                if d % k == 0:
                    v -= Fr(n[d//k], k**w)
            if v.denominator != 1:
                good = False; break
            n[d] = v.numerator
        print(f"    w={w}: " + ("INTEGRAL  " + str([n[d] for d in range(1, N)]) if good
                                else "not integral"))
