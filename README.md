# Massive Scalar Field in a Schwarzschild background

Repository to accompany our paper:

**Quasinormal modes analysis of a massive scalar field in a Schwarzschild background via the Spectral Method**

by Davide Batic, Anna Chrysostomou, Alan S. Cornell, and Denys Dutykh — *Preprint, April 2026*

- Davide Batic — Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE
- Anna Chrysostomou — Laboratoire de Physique Théorique et Hautes Énergies (LPTHE), Sorbonne Université, CNRS, Paris, France
- Alan S. Cornell — Department of Physics, University of Johannesburg, Auckland Park, South Africa
- Denys Dutykh — Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE

Here we collect the main routines used in our computations along with the raw unprocessed results.

![QNM illustration](assets/QNM_illustration.png)

## Repository structure

```
massive-scalar/
├── LICENSE
├── README.md
├── assets/
│   └── QNM_illustration.png
├── data/                          — Raw computational results (47 parameter cases)
│   └── mu_<mu>_L_<L>/            — one directory per (mu, L) pair
│       ├── raw/
│       │   ├── eigs_180.dat
│       │   ├── eigs_190.dat
│       │   └── eigs_200.dat
│       └── *_report.txt           — convergence report (when available)
├── maple/                         — Maple routines
│   ├── 02-Adapting to Spectral Method.mw
│   └── matrixassembler.mpl        — Matrix assembly using Chebyshev collocation
├── mathematica/                   — Mathematica notebooks
│   └── Massive-Scalar-QNMs_Schwarzschild_02-04-2026.nb
├── matlab/                        — MATLAB routines
│   ├── chaseeigs.m                — Polynomial eigenvalue solver for QNMs
│   └── 03-chaseeigs_sheet_resolved.m
│                                  — sheet-resolved revision (keeps Lambda and k)
├── python/                        — Independent Python cross-checks (mpmath)
│   ├── 04-massive_mode_benchmark.py — Leaver/continued-fraction benchmark
│   └── 05-potential_check.py        — high-precision effective-potential check
└── pdf/                           — Paper preprint
    └── DB-AC2-DD-QNMsMassiveScalarField-2026.pdf
```

## Computational pipeline

The routines are meant to be used in the following order; the numeric prefixes
reflect the position of each newer routine in that workflow.

### 1. Matrix assembly — `maple/matrixassembler.mpl`

`MatrixAssembler(d, n, mu, L, path)` discretises the radial equation with a
Chebyshev collocation scheme and writes the eight matrices `M0`, ..., `M7` of
the degree-7 polynomial eigenvalue problem to `assemble/M<i>_<n>.mat`, using `d`
decimal digits and `n` Chebyshev modes. The worksheet
`02-Adapting to Spectral Method.mw` documents the reduction of the massive
scalar equation to the form consumed by the assembler.

### 2. Eigenvalue extraction — `matlab/chaseeigs.m`

Reads the resolutions from `resolutions.txt` and the mass parameter from
`params.txt`, solves the polynomial eigenvalue problem with Advanpix
multiprecision `polyeig` at each resolution, maps the uniformising eigenvalue
Lambda to the frequency

```
Omega = (Lambda^2 + mu^2) / (2 Lambda),
```

discards the spurious infinite/NaN roots, and writes `results/eigs_<n>.dat` with
two columns, `Re(Omega)` and `Im(Omega)`. These are exactly the files archived
under `data/`.

### 3. Sheet-resolved eigenvalues — `matlab/03-chaseeigs_sheet_resolved.m`

The two-column output above is not sufficient for an a posteriori
classification of the modes: Lambda and `mu^2/Lambda` give the same Omega but
opposite

```
k = (Lambda^2 - mu^2) / (2 Lambda),
```

i.e. they sit on the two sheets of `k = sqrt(Omega^2 - mu^2)`. This revision
therefore keeps Lambda and k alongside Omega. It still writes the legacy
`results/eigs_<n>.dat` for backward compatibility, and additionally produces
`results/eigs_sheet_resolved_<n>.dat` with the columns

```
ReLambda ImLambda ReOmega ImOmega Rek Imk absLambda sheet_flag spatial_flag
```

where

- `sheet_flag` = `+1` on the exterior/outgoing sheet (`|Lambda| > mu`), `-1` on
  the reciprocal sheet (`|Lambda| < mu`), `0` on the branch circle;
- `spatial_flag` = `+1` for exponential spatial decay (`Im k > 0` in
  `exp(2 i k x)`), `-1` for exponential growth (`Im k < 0`), `0` when neutral.

### 4. Independent benchmark — `python/04-massive_mode_benchmark.py`

A cross-check that shares no code with the Chebyshev assembly. It evaluates the
massive-scalar three-term recurrence of R. A. Konoplya and A. V. Zhidenko,
*Phys. Lett. B* **609**, 377–384 (2005), [arXiv:gr-qc/0411059], by backward
recurrence with the Nollert large-`n` tail
`R_n = C0 + C1/sqrt(n) + C2/n`, in 80-digit `mpmath` arithmetic and in the
conventions of the paper (`M = 1`, `Omega = M omega`, `mu = M m`). The two
asymptotic branches are exposed explicitly: `out` (analytic continuation from
`k ~ +Omega` at large `|Omega|`) and `dec` (`Im k > 0`, spatially decaying
quasi-resonances). The script is a benchmark seeded near the archived spectral
roots, not a global root finder, and it reports three things:

- a direct ESM vs. CFM comparison for twelve modes with mu = 0.1, 0.2, 0.7, 0.8;
- the truncation sequence (`N_CF` = 1200, 2000, 3000, 5000) of the slowly
  converging mu = 0.2, l = 0, n = 0 root;
- a calibration against the published Hill-determinant values of Alves, Ponquio
  and Medeiros, *Phys. Rev. D* **112**, 124007 (2025).

The report is printed to the standard output and also written to a file; the
output path near the end of `main()` is hard-coded and may need to be adjusted
to a local directory before running.

### 5. Effective-potential check — `python/05-potential_check.py`

A short 60-digit `mpmath` verification of the l = 2 (i.e. L = l(l+1) = 6)
statements made about the effective potential

```
V(x) = (1 - 1/x) ( L/(4 x^2) + 1/(4 x^3) + mu^2 ).
```

It evaluates the closed-form critical mass mu_c, locates the local maximum x_M
and local minimum x_m of the potential at mu = 0.7 as roots of the cubic
`4 mu^2 x^3 - 2 L x^2 + 3 (L-1) x + 4`, compares `V(x_M)` and `V(x_m)` with
`V_infinity = mu^2`, and confirms the closed forms of the threshold point
`x_b = (5 + sqrt(43))/6` and of `mu_b^2 = -130/27 + 43 sqrt(43)/54`, at which the
local minimum touches the asymptotic value.

### 6. Figures — `mathematica/Massive-Scalar-QNMs_Schwarzschild_02-04-2026.nb`

The notebook (by A. Chrysostomou) generates Figures 1 and 2 of the paper.
Figure 1 shows the effective quasinormal mode potential as a function of
x = r/(2M) for increasing values of the scaled mass parameter mu. Figure 2 is a
two-panel figure: the left panel displays the parameter space where
V_peak > mu^2, with two benchmark points marked (mu = 0.5, L = 2 and mu = 0.5,
L = 6, where L = l(l+1)); the right panel sketches the corresponding effective
potentials.

## Requirements

- **Maple** for the matrix assembler.
- **MATLAB** with the [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/)
  for the polynomial eigenvalue solvers. Both `.m` files add the toolbox from a
  hard-coded path (`/home/dds/Soft/advanpix/`) that should be edited to match
  your installation, and expect `resolutions.txt`, `params.txt`, `assemble/` and
  `results/` in the working directory.
- **Python 3** with [`mpmath`](https://mpmath.org/) for the two cross-checks.
- **Mathematica** for the figure notebook.

## Data

The `data/` directory contains 47 parameter cases covering the following values:

- **mu** (scaled mass parameter): 0, 0.1, 0.2, 0.5, 0.7, 0.8, 1.0, 1.3, 1.4, 1.42, 1.45, 1.50, 2.0, 2.2, 3.0, 10, 100
- **L** (angular momentum, L = l(l+1)): 0, 1, 2, 10

Each case contains eigenvalues computed at 180, 190, and 200 digits of precision. Most directories also include a convergence report generated by QNM Analyser.

The archived `eigs_<n>.dat` files are in the legacy two-column format
(`Re(Omega)`, `Im(Omega)`) produced by `matlab/chaseeigs.m`; the sheet-resolved
nine-column files described in step 3 above are regenerated on demand rather
than stored here.

## Analysing the data

The raw eigenvalue files in the `data/` directory can be loaded directly into **QNM Analyser**, an interactive web dashboard for exploring convergence of quasi-normal mode eigenvalues computed at different numerical resolutions.

- Source code: <https://github.com/dutykh/qnm-analyser/>
- Live instance: <https://www.qnm-anal.denys-dutykh.com/>

For each parameter pair (mu, L), upload the three resolution files (`eigs_180.dat`, `eigs_190.dat`, `eigs_200.dat`) into QNM Analyser to automatically identify converged QNMs, classify them (general, purely imaginary, or purely real), and export publication-ready plots and reports.

## License

This project is distributed under the [GNU Lesser General Public License v2.1](LICENSE). See the `LICENSE` file for details.

## Citation

If you use the codes, routines, or data provided in this repository, please acknowledge our work by citing the following paper:

```bibtex
@article{Batic2026,
  author  = {Batic, Davide and Chrysostomou, Anna and Cornell, Alan S. and Dutykh, Denys},
  title   = {Quasinormal modes analysis of a massive scalar field in a {S}chwarzschild background via the {S}pectral {M}ethod},
  year    = {2026},
  note    = {Preprint, April 2026}
}
```
