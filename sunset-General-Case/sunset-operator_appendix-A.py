from fractions import Fraction as Fn
from math import factorial, comb

# (j,k) denotes y**j theta**k, with y on the left.  Signs are already in the
# recursion (Definition labelled def:operator), so no (-1)**j is applied at the end.
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


print(sunset_operator(4))
print(sunset_operator(5))
