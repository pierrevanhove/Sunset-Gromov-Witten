# Reproducible interval certificate (four loops)

This folder contains the computer-assisted proof of the theorem *Certified individual singularities* of [The three- and four-loop
equal-mass sunset integrals and local Calabi--Yau fourfolds and fivefolds](https://arxiv.org/pdf/2610.10821).

Let $W_j$ be the ambient solutions of $P_4W_j=\delta_{j4}$, with endpoint expansion $W_j=x(a_j+b_j\log x)+O(x^2|\log x|)$. Define

$$d_j=\frac{a_jb_0-a_0b_j}{b_0^2},\qquad
E_k=\sum_{h=1}^{k-1}\frac{(-R_\infty)^{k-1-h}}{(k-1-h)!}\,d_h .$$

The certificate proves $E_2>0$, $E_3<0$, $E_4>0$ and $E_5>0$.

## Files

| File | Role |
|---|---|
| `certify.py` | Real-ball (Arb) continuation from the germ at $y=1/512$ to the endpoint, Frobenius matching at $y=32$ or $48$, and a validated linear solve. Writes `certificate_<match>.json`. |
| `validate.py` | Exact rational checks: the 50 square-free germ coefficients, the 10 germ base cases, the induction constant $15497623/21870000$, the Frobenius bounds $1807/405$, $871/450$, $13/1350$, and the signs and overlap of the exported rational $E_k$ intervals from both certificates. Writes `validation.json` and `exact_germ_check.json`. |
| `certificate_32.json`, `certificate_48.json`, `validation.json`, `exact_germ_check.json` | Outputs of a run of the three commands below. |
| `sunset-II-certificate-worksheet.ipynb` | Worked worksheet, stored with outputs. It follows the paper section by section and checks every number the paper quotes against the program output. |
| `README.txt` | README of the Github ancillary files, with SHA-256 checksums. |

`certify.py` and `validate.py` are byte-identical to the listings printed in the paper [The three- and four-loop
equal-mass sunset integrals and local Calabi--Yau fourfolds and fivefolds](https://arxiv.org/pdf/2610.10821) and have the SHA-256 checksums listed in `README.txt`. Do not edit them.

## Running

This requires Python 3.10 or later and `python-flint==0.9.0`. Keep assertions enabled, so do not use `python3 -O`. Run exactly as in the paper:

```bash
python3 certify.py                                       # 1400 bits, N=400, match y=32, 161 steps
python3 certify.py --bits 1500 --order 420 --match 48    # independent matching check
python3 validate.py
```

A fresh run reproduces `validation.json` and `exact_germ_check.json` byte for byte. It reproduces the two certificates in every field except the wall-clock `seconds`.

## What the worksheet checks

| Paper claim | Result |
|---|---|
| 161 continuation steps (first run) | 161 |
| initial tail $e_g<5.4\times10^{-309}$ | $5.3998\times10^{-309}$ |
| largest step remainder $<1.25\times10^{-218}$ (first run) | $1.2463\times10^{-218}$ |
| radii of the $E_k$ balls $<6\times10^{-139}$ (first run) | $5.80\times10^{-139}$ |
| $Q_\infty$ enclosure, Table of $E_k$ intervals | identical in both runs |

The worksheet also verifies exactly that $Q_3=225(z+1)^2(z+2)^2$ and that $p_4=(1+y)(1+9y)(1+25y)$. It checks the conversion of $P_4$ to ordinary derivatives (`certify.ODE`) and derives the transformed equation at $y=\infty$ used in Certificate III.

## Certified values (rational enclosures, width $10^{-24}$)

| Quantity | Enclosure | Sign |
|---|---|---|
| $R_\infty$ | $-1.088825123504371170391757 \le R_\infty \le -1.088825123504371170391756$ | |
| $Q_\infty=e^{R_\infty}$ | $0.336611738699193750820017 \le Q_\infty \le 0.336611738699193750820018$ | |
| $E_2$ | $[1.959219916247366530180947,\ 1.959219916247366530180948]$ | $>0$ |
| $E_3$ | $[-0.535347537401972845058210,\ -0.535347537401972845058209]$ | $<0$ |
| $E_4$ | $[0.420379682614378171767897,\ 0.420379682614378171767898]$ | $>0$ |
| $E_5$ | $[0.748770280780437029430306,\ 0.748770280780437029430307]$ | $>0$ |

Consequently each diagonal descendant potential has an unbounded second derivative at $Q_\infty$, and the individual $Q$-Taylor radii satisfy $1/20\le\mathcal R_k\le Q_\infty$.

This is a reproducible computer-assisted proof with the analytic majorants
stated in the paper and the standard correctness assumptions for
Python/FLINT/Arb; it is not a formal proof-assistant verification.
Agreement of the two matching runs is a cross-check, not a substitute for the
error bounds.
