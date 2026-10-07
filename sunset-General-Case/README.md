# Paper I — the all-loop sunset differential operator

This folder contains the finite operator algorithm of [Equal-mass sunset integrals and genus-zero local Gromov--Witten theory of Calabi--Yau
$n$-folds](../Papers/sunset-I-general.pdf): Definition 8.1<!--`def:operator`-->, Theorem 8.2 <!--`thm:PF`-->, and Appendix A, *A symbolic implementation of the operator algorithm*. The code uses only the Python standard library and exact rational arithmetic.

| File | Role |
|---|---|
| `sunset-operator_appendix-A.py` | The code exactly as printed in Appendix A. It prints the raw coefficient dictionaries for $n=4,5$. |
| `sunset_differential_operator.py` | The same algorithm, unchanged, plus routines that print the result as a differential operator (plain text and LaTeX). It also checks the leading coefficient (8.3<!--`eq:leadP`-->) and checks that the operator annihilates the sunset period. |
| `sunset-differential-operator.ipynb` | Worked worksheet, stored with outputs. |

An operator is stored as a dictionary `(j, k) -> c_jk`, representing $\sum c_{jk}\,y^j\theta^k$ with $y$ on the left and $\theta=y\,\partial_y$. `sunset_operator(n)` returns the operator $P_{n-1}$ for the $(n-1)$-loop sunset, where $n$ is the number of propagators.

## Running

```bash
python3 sunset-operator_appendix-A.py                      # the paper's listing
python3 sunset_differential_operator.py                    # n = 4, 5
python3 sunset_differential_operator.py 2 3 4 5 6 7        # any n >= 2
python3 sunset_differential_operator.py 4 5 --latex        # also print LaTeX
python3 sunset_differential_operator.py 5 --latex --expanded
```

For $n=4$ the output is equation (8.4<!--e:P3-->) of the paper:

$$P_3=\theta^3+y(20\theta^3+30\theta^2+18\theta+4)+64y^2(\theta+1)^3 .$$

For $n=5$ (four loops) it is

$$P_4=\theta^4 + y\,Q_1(\theta) + y^2 Q_2(\theta) + y^3 Q_3(\theta),$$

$$Q_1=35\theta^4+70\theta^3+63\theta^2+28\theta+5,\quad
Q_2=259\theta^4+1036\theta^3+1580\theta^2+1088\theta+285,\quad
Q_3=225(\theta+1)^2(\theta+2)^2 .$$

This $P_4$ is the operator used by the certificate in [`../sunset-Loop-3-4`](../sunset-Loop-3-4), where it is the array `Q` in `certify.py`.

## The worksheet

`sunset-differential-operator.ipynb` covers:
1. it runs the Appendix A listing and checks that the module reproduces it for $n=2,\dots,9$;
2. it compares $n=4,5$ with Example 8.5<!--`ex:operators`-->;
3. it checks the leading coefficient $\prod_{0\le j<n/2}(1+(n-2j)^2y)$ of Theorem 8.2<!--`thm:PF`--> for $n=2,\dots,9$, including the $n=6$ value $1+56y+784y^2+2304y^3$ quoted in the text;
4. it checks that $P_{n-1}$ annihilates $\sum_D(-1)^DA_Dy^D$, a consistency check that is not part of the proof;
5. it treats $n=2,3$, which are outside the paper's standing hypothesis $n\ge4$; there the algorithm gives the one-loop and two-loop operators;
6. it displays the operators for $n=2,\dots,8$.

> The code implements a proved finite construction; it is not used to guess an
> operator from a series. The all-degree proof of the normalization is in the paper.
