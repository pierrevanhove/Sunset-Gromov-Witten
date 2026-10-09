[![Paper I: 2610.10828](https://img.shields.io/badge/Paper_I-2610.10828-b31b1b.svg)](#citation)
[![Paper II: 2610.10821](https://img.shields.io/badge/Paper_II-2610.10821-b31b1b.svg)](#citation)
[![Author: Hartmut Maennel](https://img.shields.io/badge/author-Hartmut_Maennel-blue)](#citation)
[![Author: Pierre Vanhove](https://img.shields.io/badge/author-Pierre_Vanhove-blue)](https://pierrevanhove.github.io)
![Language: Python](https://img.shields.io/badge/Language-Python-yellow?logo=python)

# Gromov–Witten theory of the multiloop sunset graph

Code accompanying two papers [Equal-mass sunset integrals and genus-zero local Gromov--Witten theory of Calabi--Yau
$n$-folds](Papers/sunset-I-general.pdf) and [The three- and four-loop
equal-mass sunset integrals and local Calabi--Yau fourfolds and fivefolds](Papers/sunset-II-loops34.pdf) by **Hartmut Maennel** and **Pierre Vanhove** relating the
equal-mass two-dimensional sunset Feynman integral to the genus-zero local
Gromov–Witten theory of the local Calabi–Yau $n$-folds $X_n=\mathrm{Tot}(K_{F_n})$,
with $F_n\subset(\mathbb P^1)^n$ a smooth $(1,\dots,1)$ hypersurface.

## Repository layout

| Path | Paper | What it does | Dependencies |
|---|---|---|---|
| [`sunset-General-Case/`](sunset-General-Case) | **[Paper I](Papers/sunset-I-general.pdf)** | The finite operator algorithm of Appendix A (*A symbolic implementation of the operator algorithm*): the all-loop sunset differential operator $P_{n-1}$, compared with Example 8.5 <!--`ex:operators`--> and Theorem 8.2<!--`thm:PF`-->. | Python standard library |
| [`sunset_GW_general_information.py`](sunset_GW_general_information.py), [`Sunset-Gromov-Witten.ipynb`](Sunset-Gromov-Witten.ipynb) | **[Paper I](Papers/sunset-I-general.pdf)**  | Exploratory companion for $L=2,\dots,6$ loops. It recovers the Picard–Fuchs operator from the sunset numbers by exact linear algebra, checks the mirror identity $f_{L-1}=(-1)^L\theta_Q F_{\gamma_L}$ exactly to $Q^{15}$, lists the invariants $N_m(\gamma_L)$, and tests integrality of the multiple-cover transforms. | Python standard library, `numpy` |
| [`sunset-Loops-3-4/`](sunset-Loops-3-4) | **[Paper II](Papers/sunset-II-loops34.pdf)** | Proof of the four-loop endpoint non-analyticity: a rigorous ball-arithmetic certificate of the signs of $E_2,\dots,E_5$ (Appendix *Reproducible interval certificate*), with a worked worksheet. | `python-flint==0.9.0` |

## Quick start

```bash
git clone https://github.com/pierrevanhove/Sunset-Gromov-Witten.git
cd Sunset-Gromov-Witten
python3 -m pip install -r requirements.txt

# Paper I: sunset operators for n = 4 (three loops) and n = 5 (four loops)
python3 sunset-General-Case/sunset_differential_operator.py 4 5 --latex

# Paper I, companion: operators, mirror identity and invariants for L = 2..6
python3 -c "from sunset_GW_general_information import report; [report(L) for L in (2,3,4,5,6)]"

# Paper II: the four-loop certificate (about 15 s), exactly as in the paper
cd sunset-Loops-3-4
python3 certify.py
python3 certify.py --bits 1500 --order 420 --match 48
python3 validate.py
```

Python 3.10 or later is required. Do **not** run with `python3 -O`, because the certificate relies on `assert` statements.

## Notebooks

All three notebooks are stored with their outputs, so they can be read on GitHub without running anything. To rerun one, start Jupyter **in the folder that contains it**, because each notebook imports the scripts next to it.

| Notebook | Content |
|---|---|
| [`sunset-General-Case/sunset-differential-operator.ipynb`](sunset-General-Case/sunset-differential-operator.ipynb) | Runs the Appendix A listing; compares with Example 8.5<!--`ex:operators`--> and with the leading coefficient (8.3<!--`eq:leadP`-->); prints the operators for $n=2,\dots,8$ in LaTeX. |
| [`Sunset-Gromov-Witten.ipynb`](Sunset-Gromov-Witten.ipynb) | `report(L)` for $L=2,\dots,6$: operator, mirror identity, $N_m(\gamma_L)$, multiple-cover integrality. |
| [`sunset-Loops-3-4/sunset-II-certificate-worksheet.ipynb`](sunset-Loops-3-4/sunset-II-certificate-worksheet.ipynb) | [Paper II](Papers/sunset-II-loops34.pdf) Certificates I–III step by step; both certificate runs; every number quoted in the paper checked against the output. |

## Conventions

* $n$ is the number of propagators, and $L=n-1$ is the number of loops.
  * [Paper I](Papers/sunset-I-general.pdf) and `sunset-General-Case/` use $n$; `sunset_GW_general_information.py` uses $L$.
  * The Fano variety is therefore $F_n\subset(\mathbb P^1)^n$ in the papers and $F_L\subset(\mathbb P^1)^{L+1}$ in the script.
* $y=-t$ throughout [Paper I](Papers/sunset-I-general.pdf), Remark 8.3<!--`rem:unsigned`-->). The holomorphic period is $\sum_D(-1)^DA_Dy^D$, where $A_D=\sum_{|d|=D}\binom{D}{d_1,\dots,d_n}^2$.
* Operators are written $\sum_j y^jP_j(\theta_y)$ with $\theta_y=y\,\partial_y$ and powers of $y$ on the left.

## Citation

```bibtex
@article{Maennel:2026lud,
    author = "Maennel, Hartmut and Vanhove, Pierre",
    title = "{Equal-mass sunset integrals and genus-zero local Gromov--Witten theory of Calabi--Yau $n$-folds}",
    eprint = "2610.10828",
    archivePrefix = "arXiv",
    primaryClass = "math.AG",
    month = "10",
    year = "2026"
}
@article{Maennel:2026fhs,
    author = "Maennel, Hartmut and Vanhove, Pierre",
    title = "{The three- and four-loop equal-mass sunset integrals and local Calabi--Yau fourfolds and fivefolds}",
    eprint = "2610.10821",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    month = "10",
    year = "2026"
}
```

See also P. Vanhove, *The multiloop sunset to all orders*, Ann. Henri Poincaré (2026),
[arXiv:2603.03183](https://arxiv.org/abs/2603.03183).
