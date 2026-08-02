#!/usr/bin/env python3
"""Reproduce the massive-mode benchmarks added in referee issue 7.

The independent recurrence is the massive-scalar three-term relation of
R. A. Konoplya and A. V. Zhidenko, Phys. Lett. B 609, 377-384 (2005),
arXiv:gr-qc/0411059.  We work in the manuscript convention M=1,
Omega=M*omega, mu=M*m, and evaluate the continued fraction by backward
recurrence with the Nollert large-n tail quoted in that paper.

Two asymptotic branches are exposed explicitly:

  out : analytic continuation from k(Omega) ~ +Omega at large |Omega|;
  dec : Im k > 0 in exp(2 i k x), giving spatial decay as x -> infinity.

The code is intentionally independent of the Chebyshev matrix assembly.
It is a benchmark seeded near the archived spectral roots, not a global
root finder.

Authors: Davide Batic (Khalifa University of Science and Technology,
         Abu Dhabi, UAE)
         Anna Chrysostomou (LPTHE, Sorbonne Universite, CNRS, Paris, France)
         Alan S. Cornell (University of Johannesburg, Auckland Park,
         South Africa)
         Dr. Denys Dutykh (Khalifa University of Science and Technology,
         Abu Dhabi, UAE)
Last modified: 2 August 2026
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

import mpmath as mp

mp.mp.dps = 80
Branch = Literal["out", "dec"]


@dataclass(frozen=True)
class Case:
    mu: mp.mpf
    ell: int
    index: str
    branch: Branch
    guess: mp.mpc
    omega_spectral: mp.mpc
    truncations: tuple[int, ...] = (300, 500, 800, 1200)


def z(re: str, im: str) -> mp.mpc:
    return mp.mpc(mp.mpf(re), mp.mpf(im))


def choose_k(omega: mp.mpc, mu: mp.mpf, branch: Branch) -> mp.mpc:
    """Choose k=sqrt(Omega^2-mu^2) on the requested asymptotic branch."""
    k = mp.sqrt(omega * omega - mu * mu)
    if branch == "out":
        # Konoplya-Zhidenko continuation: k occupies the same quadrant as
        # Omega, implementing k~+Omega away from the branch points.
        if mp.re(omega) >= 0 and mp.re(k) < 0:
            k = -k
        if mp.re(omega) < 0 and mp.re(k) > 0:
            k = -k
        if mp.im(omega) < 0 and mp.im(k) > 0:
            k = -k
        if mp.im(omega) > 0 and mp.im(k) < 0:
            k = -k
    elif branch == "dec":
        # In the manuscript asymptotic factor exp(2 i k x), Im k>0
        # produces exponential spatial decay.
        if mp.im(k) < 0:
            k = -k
    else:  # pragma: no cover
        raise ValueError(f"unknown branch {branch!r}")
    return k


def recurrence_coefficients(
    n: int, omega: mp.mpc, mu: mp.mpf, ell: int, k: mp.mpc
) -> tuple[mp.mpc, mp.mpc, mp.mpc]:
    """Konoplya-Zhidenko three-term recurrence in units M=1."""
    wp = omega + k
    alpha = (n + 1) * (n + 1 - 4j * omega)
    beta = (
        wp * (4 * wp**2 + 1j * (2 * n + 1) * (omega + 3 * k)) / k
        - 2 * n * (n + 1)
        - 1
        - ell * (ell + 1)
    )
    gamma = (n - 1j * wp**2 / k) ** 2
    return alpha, beta, gamma


def nollert_tail(n: int, k: mp.mpc, mu: mp.mpf) -> mp.mpc:
    """Nollert remainder R_n=C0+C1/sqrt(n)+C2/n through O(n^-1)."""
    c0 = mp.mpf(-1)
    c1 = 2 * mp.sqrt(-1j * k)
    if mp.re(c1) < 0:
        c1 = -c1
    c2 = mp.mpf(3) / 4 + 4j * k + 1j * mu**2 / k
    nn = mp.mpf(n)
    return c0 + c1 / mp.sqrt(nn) + c2 / nn


def continued_fraction(
    omega: mp.mpc, mu: mp.mpf, ell: int, truncation: int, branch: Branch
) -> mp.mpc:
    k = choose_k(omega, mu, branch)
    remainder = nollert_tail(truncation + 1, k, mu)
    for n in range(truncation, 0, -1):
        alpha, beta, gamma = recurrence_coefficients(n, omega, mu, ell, k)
        remainder = gamma / (beta - alpha * remainder)
    alpha0, beta0, _ = recurrence_coefficients(0, omega, mu, ell, k)
    return beta0 - alpha0 * remainder


def solve(case: Case, truncation: int, seed: mp.mpc | None = None) -> mp.mpc:
    f = lambda omega: continued_fraction(
        omega, case.mu, case.ell, truncation, case.branch
    )
    omega0 = case.guess if seed is None else seed
    omega1 = omega0 + z("1e-7", "-1e-7")
    return mp.findroot(f, (omega0, omega1), tol=mp.mpf("1e-45"), maxsteps=100)


def solve_sequence(case: Case) -> tuple[mp.mpc, mp.mpf, list[mp.mpc]]:
    roots: list[mp.mpc] = []
    seed = case.guess
    for truncation in case.truncations:
        root = solve(case, truncation, seed=seed)
        roots.append(root)
        seed = root
    drift = max(abs(a - b) for a in roots for b in roots)
    return roots[-1], drift, roots


def direct_cases() -> Sequence[Case]:
    return [
        Case(mp.mpf("0.1"), 0, "n=0", "out", z("0.1123610565", "-0.0968230289"),
             z("0.11236105649688380", "-0.096823028935690725")),
        Case(mp.mpf("0.1"), 1, "n=0", "out", z("0.2974156612", "-0.0949570736"),
             z("0.29741566124545843", "-0.094957073606208892")),
        Case(mp.mpf("0.1"), 1, "n=1", "out", z("0.2646885510", "-0.3028507124"),
             z("0.26468855096741711", "-0.30285071241979511")),
        Case(mp.mpf("0.1"), 2, "n=0", "out", z("0.4868037520", "-0.0956745815"),
             z("0.48680375196030301", "-0.095674581470046954")),
        Case(mp.mpf("0.1"), 2, "n=2", "out", z("0.4306647741", "-0.5064640671"),
             z("0.43066477411230231", "-0.50646406706345526")),
        Case(mp.mpf("0.2"), 0, "n=0", "out", z("0.11620657", "-0.07535829"),
             z("0.11620668013418312", "-0.075358181321549700"),
             (1200, 2000, 3000, 5000)),
        Case(mp.mpf("0.2"), 1, "n=0", "out", z("0.3109569084", "-0.0865932856"),
             z("0.31095690839826745", "-0.086593285615637949")),
        Case(mp.mpf("0.2"), 2, "n=0", "out", z("0.4963266059", "-0.0923891664"),
             z("0.49632660594248351", "-0.092389166442148257")),
        Case(mp.mpf("0.7"), 0, "j=0", "dec", z("0.6954459954", "-0.0019453419"),
             z("0.69544599544276420", "-0.0019453418680002625")),
        Case(mp.mpf("0.7"), 0, "j=1", "dec", z("0.6937304671", "-0.0032365295"),
             z("0.69373046713287423", "-0.0032365294751243311")),
        Case(mp.mpf("0.8"), 0, "j=0", "dec", z("0.7935418414", "-0.0032252458"),
             z("0.79354184141025541", "-0.0032252458168555735")),
        Case(mp.mpf("0.8"), 1, "j=0", "dec", z("0.7822884986", "-0.0101622674"),
             z("0.78228849861756211", "-0.010162267390471933")),
    ]


def format_complex(value: mp.mpc, digits: int = 26) -> str:
    re = mp.nstr(mp.re(value), digits)
    im = mp.nstr(abs(mp.im(value)), digits)
    sign = "+" if mp.im(value) >= 0 else "-"
    return f"{re} {sign} {im} i"


def main() -> None:
    lines: list[str] = []
    lines.append("Independent massive-scalar Leaver-Nollert benchmark")
    lines.append("Arithmetic: 80 decimal digits")
    lines.append("")
    lines.append("DIRECT ESM-CFM COMPARISON")
    lines.append("mu ell index branch Omega_ESM Omega_CFM |Delta| CF_max_drift CF_final_step")

    for case in direct_cases():
        root, drift, roots = solve_sequence(case)
        final_step = abs(roots[-1] - roots[-2])
        discrepancy = abs(case.omega_spectral - root)
        lines.append(
            f"{mp.nstr(case.mu, 4):>4} {case.ell:>3} {case.index:>5} "
            f"{case.branch:>4}  ESM=({format_complex(case.omega_spectral)})  "
            f"CFM=({format_complex(root)})  "
            f"Delta={mp.nstr(discrepancy, 12)}  "
            f"CF-max-drift={mp.nstr(drift, 12)}  "
            f"CF-final-step={mp.nstr(final_step, 12)}  trunc={case.truncations}"
        )

    lines.append("")
    lines.append("SLOW-ROOT TRUNCATION SEQUENCE: mu=0.2, ell=0, n=0, out")
    slow_case = next(
        case for case in direct_cases()
        if case.mu == mp.mpf("0.2") and case.ell == 0 and case.index == "n=0"
    )
    _, _, slow_roots = solve_sequence(slow_case)
    for truncation, root in zip(slow_case.truncations, slow_roots):
        lines.append(f"N_CF={truncation}: Omega=({format_complex(root, 32)})")
    lines.append(
        "N_CF=3000 -> 5000 change="
        + mp.nstr(abs(slow_roots[-1] - slow_roots[-2]), 16)
    )

    lines.append("")
    lines.append("CALIBRATION AGAINST PUBLISHED HILL-DETERMINANT TABLE")
    lines.append("Ref.: Alves, Ponquio, Medeiros, Phys. Rev. D 112, 124007 (2025)")
    lines.append("Their printed m=0.1 maps to mu=0.05 after matching the radial mass term.")
    published = {
        0: z("0.110988", "-0.102842"),
        1: z("0.294054", "-0.096988"),
        2: z("0.484433", "-0.0964882"),
    }
    guesses = {
        0: z("0.110995", "-0.102846"),
        1: z("0.294054", "-0.096988"),
        2: z("0.484433", "-0.096488"),
    }
    for ell in (0, 1, 2):
        case = Case(mp.mpf("0.05"), ell, "n=0", "out", guesses[ell], z("0", "0"))
        root, drift, roots = solve_sequence(case)
        final_step = abs(roots[-1] - roots[-2])
        difference = abs(root - published[ell])
        lines.append(
            f"ell={ell}: CFM=({format_complex(root)})  "
            f"published=({format_complex(published[ell], 16)})  "
            f"difference={mp.nstr(difference, 12)}  "
            f"CF-max-drift={mp.nstr(drift, 12)}  "
            f"CF-final-step={mp.nstr(final_step, 12)}"
        )

    output = "\n".join(lines) + "\n"
    print(output, end="")
    Path("/mnt/data/DR13837_issue7_benchmark_results.txt").write_text(
        output, encoding="utf-8"
    )


if __name__ == "__main__":
    main()
