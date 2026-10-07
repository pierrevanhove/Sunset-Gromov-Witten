"""
Finite operator algorithm of Definition def:operator (sunset-I-general.tex).

    n   = number of propagators (paper notation);  L = n-1 loops (sunset.py notation)
    P_{n-1}(y, theta) = theta^{n-1} + sum_{j>=1} y^j P_j(theta),   theta = y d/dy,
    powers of y on the left, y = -t (Remark rem:unsigned).

The first block (add ... sunset_operator) is the code listed in Appendix app:code,
unchanged.  The second block turns the coefficient dictionary into the differential
operator, in the same format as sunset.py, and in LaTeX as in Example ex:operators.

The algorithm also runs for n = 2, 3 (one and two loops), which lie outside the
standing hypothesis n >= 4 of the paper.

Usage:
    python sunset_differential_operator.py                   # n = 4, 5
    python sunset_differential_operator.py 2 3 4 5 6 7       # any n >= 2
    python sunset_differential_operator.py 4 5 --latex       # also print LaTeX
    python sunset_differential_operator.py 5 --latex --expanded
"""
from fractions import Fraction as Fn
from math import factorial, comb, gcd

# ======================================================================
# Appendix app:code -- the finite operator algorithm (unchanged)
# (j,k) denotes y**j theta**k, with y on the left.  Signs are already in the
# recursion (Definition labelled def:operator), so no (-1)**j is applied at the end.
# ======================================================================
def add(a,b):
    c=a.copy()
    for key,v in b.items(): c[key]=c.get(key,Fn(0))+v
    return {key:v for key,v in c.items() if v}
def scale(a,c): return {key:c*v for key,v in a.items() if c*v}
def theta_left(a):
    c={}
    for (j,k),v in a.items():
        c=add(c,{(j,k+1):v,(j,k):j*v})
    return c
def y_left(a): return {(j+1,k):v for (j,k),v in a.items()}
def right_theta_plus(a,r):
    c={}
    for (j,k),v in a.items():
        c=add(c,{(j,k+1):v,(j,k):r*v})
    return c
def sunset_operator(n):
    T=[{(0,0):Fn(1)},{(0,1):Fn(1,n)}]
    for k in range(1,n):
        T.append(scale(add(theta_left(T[k]),
            scale(y_left(T[k-1]),k)),Fn(1,n-k)))
    U=scale(add(theta_left(T[n]),scale(y_left(T[n-1]),n)),factorial(n))
    # Exact normalization sanity check; the all-degree proof is in the text.
    assert {key:v for key,v in U.items() if key[0]==0}=={(0,n+1):Fn(1)}
    P={(0,n-1):Fn(1)}
    for j in sorted({j for j,k in U if j}):
        term={(j,k):v for (jj,k),v in U.items() if jj==j}
        for r in range(1,j):
            term=right_theta_plus(right_theta_plus(term,r),r)
        P=add(P,term)
    return P

# ======================================================================
# From coefficients to the differential operator
# ======================================================================
def operator_table(n):
    """(J, P) with P[j][k] = coefficient of y^j theta^k, k = 0..n-1.
       Same format as find_operator(L, a) in sunset.py, with L = n-1."""
    D = sunset_operator(n)
    J = max(j for j, _ in D)
    P = [[D.get((j, k), Fn(0)) for k in range(n)] for j in range(J + 1)]
    return J, [[int(c) if c.denominator == 1 else c for c in row] for row in P]

def polyfmt(row):
    """identical to polyfmt in sunset.py"""
    t = []
    for k in range(len(row) - 1, -1, -1):
        if row[k]:
            t.append(f"{row[k]}" + ("" if k == 0 else ("*th" if k == 1 else f"*th^{k}")))
    return (" + ".join(t) or "0").replace("+ -", "- ")

def leading_symbol(n):
    """coefficient of theta^{n-1} as a polynomial in y, and the prediction (eq:leadP)
       prod_{0<=j<n/2} (1 + (n-2j)^2 y)."""
    J, P = operator_table(n)
    lead = [P[j][n - 1] for j in range(J + 1)]
    pred = [1]
    for j in range((n + 1) // 2):
        m = (n - 2 * j) ** 2
        pred = [(pred[i] if i < len(pred) else 0) + (m * pred[i - 1] if i else 0)
                for i in range(len(pred) + 1)]
    return lead, pred

def sunset_numbers(n, Dmax):
    """A_D = sum_{d_1+..+d_n=D} (D!/(d_1!...d_n!))^2, by convolution."""
    g = [Fn(1, factorial(d) ** 2) for d in range(Dmax + 1)]
    s = [Fn(1)] + [Fn(0)] * Dmax
    for _ in range(n):
        s = [sum(s[i] * g[D - i] for i in range(D + 1)) for D in range(Dmax + 1)]
    return [int(factorial(D) ** 2 * s[D]) for D in range(Dmax + 1)]

def check_recurrence(n, Dmax=40):
    """P_{n-1} annihilates the holomorphic period sum_D (-1)^D A_D y^D
       (the convention of sunset.py, i.e. y = -t)."""
    J, P = operator_table(n)
    a = [(-1) ** D * A for D, A in enumerate(sunset_numbers(n, Dmax))]
    return all(sum(sum(c * (D - j) ** k for k, c in enumerate(P[j])) * a[D - j]
                   for j in range(J + 1) if D - j >= 0) == 0
               for D in range(Dmax + 1))

def report_operator(n):
    """print the operator in the same format as report() of sunset.py."""
    L = n - 1
    J, P = operator_table(n)
    lead, pred = leading_symbol(n)
    print("=" * 76)
    print(f"  n = {n} propagators,  L = {L} loops")
    print("=" * 76)
    print(f"  sunset numbers A_D : {sunset_numbers(n, 6)}")
    print(f"  Picard-Fuchs operator: J = {J}   ({J+1}-term recurrence)")
    for j in range(J + 1):
        print(f"     P_{j}(th) = {polyfmt(P[j])}")
    print(f"  coefficient of th^{L}: {lead}   (eq:leadP) {pred}:  {lead == pred}")
    print(f"  finite singular points s = -1/y : {sorted((n - 2*j)**2 for j in range((n + 1)//2))}")
    print(f"  annihilates sum (-1)^D A_D y^D to y^40 : {check_recurrence(n)}")

# ---------------------------------------------------------------- LaTeX
def _divide_linear(p, r):
    """divide p(theta) (ascending coefficients) by (theta + r); None if not exact."""
    q = [0] * (len(p) - 1); rem = p[-1]
    for k in range(len(p) - 2, -1, -1):
        q[k] = rem
        rem = p[k] - r * rem
    return q if rem == 0 else None

def _latex_poly(p, th):
    terms = []
    for k in range(len(p) - 1, -1, -1):
        c = p[k]
        if not c: continue
        mon = "" if k == 0 else (th if k == 1 else (f"{th}^{k}" if k < 10 else f"{th}^{{{k}}}"))
        cs = str(abs(c)) if (abs(c) != 1 or k == 0) else ""
        terms.append(("-" if c < 0 else "+") + cs + mon)
    s = "".join(terms)
    return s[1:] if s.startswith("+") else s

def _latex_factored(p, th, J):
    """(c, "(theta+1)^a (theta+2)^b ...") if p splits into such factors, else None."""
    while len(p) > 1 and p[-1] == 0: p = p[:-1]
    c = 0
    for x in p: c = gcd(c, x)
    p = [x // c for x in p]
    mult = {}
    for r in range(1, J + 1):
        while len(p) > 1:
            q = _divide_linear(p, r)
            if q is None: break
            p = q; mult[r] = mult.get(r, 0) + 1
    if len(p) != 1 or not mult: return None
    c *= p[0]
    fac = "".join(f"({th}+{r})" + (f"^{m}" if m > 1 else "") for r, m in sorted(mult.items()))
    return ("" if c == 1 else str(c)), fac

def latex_operator(n, factor=True):
    """P_{n-1} in the notation of Example ex:operators."""
    th = r"\theta_y"
    J, P = operator_table(n)
    out = [th if n == 2 else (f"{th}^{n-1}" if n < 11 else f"{th}^{{{n-1}}}")]
    for j in range(1, J + 1):
        yj = "y" if j == 1 else f"y^{{{j}}}"
        f = _latex_factored(P[j], th, J) if factor else None
        out.append(f"{f[0]}{yj}{f[1]}" if f else f"{yj}({_latex_poly(P[j], th)})")
    return f"P_{{{n-1}}}=" + "+".join(out)

# ======================================================================
if __name__ == "__main__":
    import sys
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ns = [int(a) for a in args] or [4, 5]
    for n in ns:
        report_operator(n)
        if "--latex" in sys.argv:
            print("\n  LaTeX:\n  " + latex_operator(n, factor="--expanded" not in sys.argv))
        print()