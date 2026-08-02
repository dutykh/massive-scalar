#!/usr/bin/env python3
"""High-precision check of the l=2 massive-scalar potential around mu=0.7.

Authors: Davide Batic (Khalifa University of Science and Technology,
         Abu Dhabi, UAE)
         Anna Chrysostomou (LPTHE, Sorbonne Universite, CNRS, Paris, France)
         Alan S. Cornell (University of Johannesburg, Auckland Park,
         South Africa)
         Dr. Denys Dutykh (Khalifa University of Science and Technology,
         Abu Dhabi, UAE)
Last modified: 2 August 2026
"""

import mpmath as mp

mp.mp.dps = 60
L = mp.mpf(6)
mu = mp.mpf('0.7')


def potential(x: mp.mpf, mass: mp.mpf) -> mp.mpf:
    return (1 - 1 / x) * (L / (4 * x**2) + 1 / (4 * x**3) + mass**2)


def derivative_polynomial(x: mp.mpf, mass: mp.mpf) -> mp.mpf:
    return 4 * mass**2 * x**3 - 2 * L * x**2 + 3 * (L - 1) * x + 4

mu_c_sq = (
    mp.sqrt(3) * (3 * L**2 + 2 * L + 3) ** (mp.mpf(3) / 2)
    - 9 * (L - 1) * (L + 1) ** 2
) / 288
mu_c = mp.sqrt(mu_c_sq)

roots = mp.polyroots([4 * mu**2, -2 * L, 3 * (L - 1), 4], maxsteps=1000)
positive_roots = sorted([mp.re(r) for r in roots if abs(mp.im(r)) < mp.mpf('1e-45') and mp.re(r) > 1])
x_max, x_min = positive_roots
v_max = potential(x_max, mu)
v_min = potential(x_min, mu)

x_b = (5 + mp.sqrt(43)) / 6
mu_b_sq = -mp.mpf(130) / 27 + 43 * mp.sqrt(43) / 54
mu_b = mp.sqrt(mu_b_sq)

print(f"mu_c^2 = {mp.nstr(mu_c_sq, 45)}")
print(f"mu_c   = {mp.nstr(mu_c, 45)}")
print(f"x_M(mu=0.7) = {mp.nstr(x_max, 45)}")
print(f"x_m(mu=0.7) = {mp.nstr(x_min, 45)}")
print(f"V_M(mu=0.7) = {mp.nstr(v_max, 45)}")
print(f"V_m(mu=0.7) = {mp.nstr(v_min, 45)}")
print(f"V_infinity   = {mp.nstr(mu**2, 45)}")
print(f"P(x_M)       = {mp.nstr(derivative_polynomial(x_max, mu), 8)}")
print(f"P(x_m)       = {mp.nstr(derivative_polynomial(x_min, mu), 8)}")
print(f"x_b          = {mp.nstr(x_b, 45)}")
print(f"mu_b^2       = {mp.nstr(mu_b_sq, 45)}")
print(f"mu_b         = {mp.nstr(mu_b, 45)}")
print(f"V(x_b,mu_b)-mu_b^2 = {mp.nstr(potential(x_b, mu_b)-mu_b_sq, 8)}")
print(f"P(x_b,mu_b)         = {mp.nstr(derivative_polynomial(x_b, mu_b), 8)}")
