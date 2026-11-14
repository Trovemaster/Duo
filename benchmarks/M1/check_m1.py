#!/usr/bin/env python3
"""Analytic and independent checks of Duo magnetic-dipole (M1) line strengths.

    python check_m1.py /path/to/duo.x

Each input in this folder runs in a temporary directory. Line strengths S are
in mu_B^2 (Duo prints them under its generic 'Debye**2' header).

A  M1_2Sigma_caseb.inp: one 2Sigma+ state with constant spin-rotation, Hund's
   case (b). g_s S cannot change N, so lines with Delta N != 0 or Delta v != 0
   vanish (e.g. J=0.5,N=0 -> J=1.5,N=2). The F1-F2 lines inside one N have
   S = g_s^2 (2J+1)(2J'+1) {S J' N; J S 1}^2 S(S+1)(2S+1).
   This fixes the sign of the spin ladders relative to Lambda + g_s Sigma.
B  M1_p_electron_cased.inp: 2Sigma+ and 2Pi of one p electron, case (d).
   L + g_s S commutes with R^2, so every line with nonzero frequency vanishes,
   in Duo and in the separate contraction below. This fixes the orbital sign.
C  M1_2Pi_2Sigma.inp: Duo strengths of all lines agree with the separate
   contraction of the eigenvectors that Duo saves (CHECKPOINT, ascii).

Separate contraction: mu = L + g_s S (the common electron sign drops out),
T(0) = Lambda + g_s Sigma, T(+-1) = -+(L+- + g_s S+-)/sqrt(2), and
S = (2J'+1)(2J+1) |sum c' c (-1)^(J'-Omega') (J' 1 J; -Omega' q Omega) <'|T_q|>|^2.
A Duo L+ field value is <higher Lambda|L+|lower Lambda> (as in the Duo L.S
and J.L Hamiltonian terms); its Lambda<0 image has (-1)^(number of Sigma-).
Tolerances: nulls |S| <= 1e-10 mu_B^2, about 1e-11 of the allowed strengths
(~1-10 mu_B^2); nonzero S agree to 1e-7 relative + 1e-10 absolute, above the
9 digits Duo prints and the 12 digits of the saved eigenvectors.
Exit status 0: PASS, 1: FAIL.
"""
import re
import subprocess
import sys
import tempfile
from math import factorial, sqrt
from pathlib import Path

GS = 2.00231930436256        # g_s used by Duo (accuracy.f90)
NULL, RTOL, ATOL = 1e-10, 1e-7, 1e-10
HERE = Path(__file__).resolve().parent
# L+ field values and Sigma- states; they must match the input files.
MODELS = {'M1_p_electron_cased': ({(2, 1): sqrt(2)}, set()),
          'M1_2Pi_2Sigma': ({(2, 1): 0.8, (2, 3): 0.6}, {3})}


def fac(x):
    return factorial(round(x))


def threej(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol, Racah formula."""
    if abs(m1 + m2 + m3) > 1e-8 or not abs(j1 - j2) <= j3 <= j1 + j2 \
            or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    kmin = round(max(0, j2 - j3 - m1, j1 - j3 + m2))
    kmax = round(min(j1 + j2 - j3, j1 - m1, j2 + m2))
    total = sum((-1)**k / (fac(k) * fac(j3 - j2 + k + m1) * fac(j3 - j1 + k - m2) * fac(j1 + j2 - j3 - k)
                           * fac(j1 - k - m1) * fac(j2 - k + m2)) for k in range(kmin, kmax + 1))
    return (-1)**round(j1 - j2 - m3) * total * sqrt(
        fac(j1 + j2 - j3) * fac(j1 - j2 + j3) * fac(-j1 + j2 + j3) / fac(j1 + j2 + j3 + 1)
        * fac(j1 + m1) * fac(j1 - m1) * fac(j2 + m2) * fac(j2 - m2) * fac(j3 + m3) * fac(j3 - m3))


def sixj(a, b, c, d, e, f):
    """Wigner 6j symbol {a b c; d e f}, Racah formula (triangles assumed)."""
    def delta(x, y, z):
        return sqrt(fac(x + y - z) * fac(x - y + z) * fac(-x + y + z) / fac(x + y + z + 1))
    tmin = round(max(a + b + c, a + e + f, d + b + f, d + e + c))
    tmax = round(min(a + b + d + e, b + c + e + f, c + a + f + d))
    total = sum((-1)**t * fac(t + 1) / (fac(t - a - b - c) * fac(t - a - e - f) * fac(t - d - b - f)
                                         * fac(t - d - e - c) * fac(a + b + d + e - t)
                                         * fac(b + c + e + f - t) * fac(c + a + f + d - t))
                for t in range(tmin, tmax + 1))
    return delta(a, b, c) * delta(a, e, f) * delta(d, b, f) * delta(d, e, c) * total


def run_duo(duo, name, folder):
    out = folder / f'{name}.out'
    with open(HERE / f'{name}.inp') as stdin, open(out, 'w') as stdout:
        subprocess.run([duo], stdin=stdin, stdout=stdout, cwd=folder, check=True)
    text = out.read_text()
    zpe = float(re.findall(r'Zero point energy \(ZPE\) =\s*(\S+)', text)[-1])
    # absorption lines: J', sym', J, sym, branch, E', E, nu, S, A, I, then the same number of
    # assignment fields for ' and for the lower level; '<-' and brackets are dropped
    lines = [[x for x in row.replace('(', ' ').replace(')', ' ').split() if x != '<-']
             for row in text.splitlines() if '<-' in row and ("A'" in row or 'A"' in row)]
    return zpe, [dict(J=(float(r[0]), float(r[2])), sym=(r[1], r[3]), E=(float(r[5]), float(r[6])),
                      S=float(r[8]), v=(int(r[12]), int(r[12 + (len(r) - 11) // 2])))
                 for r in lines]


def read_levels(folder, name, zpe):
    """Levels (J, sym, E - ZPE, coefficients), basis rows per J, radial functions."""
    rows = iter((folder / f'{name}_vib.chk').read_text().splitlines())
    next(rows)
    npoints = int(next(rows).split('=')[1])
    radial = {}
    for head in rows:
        if head.startswith('End'):
            break
        state, v = map(int, head.split()[2:4])
        radial[state, v] = [float(next(rows)) for _ in range(npoints)]
    energy = {}
    for row in (folder / f'{name}_values.chk').read_text().splitlines():
        cols = row.split()
        if len(cols) > 10 and cols[0].isdigit():
            energy[float(cols[2]), int(cols[3]), int(cols[0])] = float(cols[1]) - zpe
    basis, levels = {}, []
    rows = iter((folder / f'{name}_vectors.chk').read_text().splitlines())
    for row in rows:
        if row.startswith('basis set quantum numbers/start'):
            next(rows)
            block = []
            for row in rows:
                if row.startswith('basis set'):
                    break
                state, v, lam, spin, sigma, omega = row.split()
                block.append((int(state), int(v), int(lam), float(spin), float(sigma), float(omega)))
        elif row.startswith('J ='):
            J, sym = float(row.split()[2]), int(row.split()[5])
            basis[J] = block
        elif row.startswith('Root ='):
            root, size = int(row.split()[2]), int(row.split()[5])
            levels.append(dict(J=J, sym=["A'", 'A"'][sym - 1], E=energy[J, sym, root],
                               c=[float(next(rows)) for _ in range(size)]))
    return levels, basis, radial


def m1_matrix(bra, ket, Jf, Ji, radial, lplus, sigma_minus):
    """<f|...|i> = (-1)^(Jf-Omega_f) (Jf 1 Ji; -Omega_f q Omega_i) <f|T_q|i> on the basis rows."""
    matrix = [[0.0] * len(ket) for _ in bra]
    for f, (sf, vf, lf, spin_f, mf, of) in enumerate(bra):
        for i, (si, vi, li, spin_i, mi, oi) in enumerate(ket):
            q = round(of - oi)
            if abs(q) > 1 or spin_f != spin_i:
                continue
            if (sf, vf, lf, mf) == (si, vi, li, mi):
                t = li + GS * mi
            elif (sf, vf, lf) == (si, vi, li) and mf == mi + q:
                t = -q / sqrt(2) * GS * sqrt(spin_i * (spin_i + 1) - mi * (mi + q))
            elif mf == mi and lf == li + q and ((sf, si) in lplus or (si, sf) in lplus):
                value = lplus.get((sf, si), lplus.get((si, sf)))
                if lf + li < 0:
                    value *= (-1)**(len({sf, si} & sigma_minus))
                overlap = sum(a * b for a, b in zip(radial[sf, vf], radial[si, vi]))
                t = -q / sqrt(2) * value * overlap
            else:
                continue
            matrix[f][i] = (-1)**round(Jf - of) * threej(Jf, 1, Ji, -of, q, oi) * t
    return matrix


def contraction(levels, basis, radial, lplus, sigma_minus):
    """S for every same-symmetry pair with Delta J <= 1 and nu > 0.001 cm-1, key (upper, lower)."""
    strengths, matrices = {}, {}
    for up in levels:
        for low in levels:
            if up['sym'] != low['sym'] or abs(up['J'] - low['J']) > 1 or up['E'] - low['E'] <= 1e-3:
                continue
            key = up['J'], low['J']
            if key not in matrices:
                matrices[key] = m1_matrix(basis[key[0]], basis[key[1]], *key, radial, lplus, sigma_minus)
            amplitude = sum(cu * sum(m * cl for m, cl in zip(row, low['c']))
                            for cu, row in zip(up['c'], matrices[key]))
            strengths[id(up), id(low)] = (2 * up['J'] + 1) * (2 * low['J'] + 1) * amplitude**2
    return strengths


def check_caseb(duo_lines, *_):
    def n_rot(J, sym):  # 2Sigma+: parity (A' = +) is (-1)^N
        return J - 0.5 if (sym == "A'") == (round(J - 0.5) % 2 == 0) else J + 0.5
    spin_flip, worst_null, errors = 0, 0.0, []
    for line in duo_lines:
        (Jf, Ji), (sf, si) = line['J'], line['sym']
        N = n_rot(Jf, sf)
        if N == n_rot(Ji, si) and line['v'][0] == line['v'][1]:
            spin_flip += 1
            exact = GS**2 * (2 * Jf + 1) * (2 * Ji + 1) * sixj(0.5, Jf, N, Ji, 0.5, 1)**2 * 1.5
            if abs(line['S'] - exact) > RTOL * exact:
                errors.append(f"F1-F2 N={N:g} J={Ji}->{Jf}: S={line['S']:.8e}, analytic {exact:.8e}")
        else:
            worst_null = max(worst_null, abs(line['S']))
    if worst_null > NULL:
        errors.append(f'largest Delta N/Delta v line S = {worst_null:.3e} > {NULL:g}')
    if spin_flip != 12:  # N = 1, 2, 3 (J <= 3.5) for v = 0..3
        errors.append(f'{spin_flip} F1-F2 lines found, 12 expected')
    return errors, f'{spin_flip} F1-F2 lines match the 6j strengths, other lines |S| <= {worst_null:.1e}'


def check_cased(duo_lines, levels, strengths):
    worst_duo = max((abs(line['S']) for line in duo_lines), default=0.0)
    worst_own = max(strengths.values(), default=0.0)
    errors = [f'{who} largest S = {s:.3e} > {NULL:g}'
              for who, s in (('Duo', worst_duo), ('separate contraction', worst_own)) if s > NULL]
    return errors, (f'{len(duo_lines)} Duo lines |S| <= {worst_duo:.1e}, '
                    f'{len(strengths)} contracted pairs |S| <= {worst_own:.1e}')


def check_compare(duo_lines, levels, strengths):
    def level(J, sym, E):
        found = [x for x in levels if x['J'] == J and x['sym'] == sym and abs(x['E'] - E) < 2e-4]
        if len(found) != 1:
            raise ValueError(f'{len(found)} levels match J={J} {sym} E={E}')
        return found[0]
    errors, matched, worst = [], set(), 0.0
    for line in duo_lines:
        key = tuple(id(level(J, s, E)) for J, s, E in zip(line['J'], line['sym'], line['E']))
        matched.add(key)
        own = strengths.get(key, 0.0)
        worst = max(worst, abs(line['S'] - own) / max(own, ATOL))
        if abs(line['S'] - own) > ATOL + RTOL * own:
            errors.append(f"J={line['J'][1]}->{line['J'][0]} E={line['E']}: Duo {line['S']:.8e}, separate {own:.8e}")
    missing = [s for key, s in strengths.items() if key not in matched and s > ATOL]
    if missing:
        errors.append(f'{len(missing)} contracted lines with S > {ATOL:g} are missing in Duo, largest {max(missing):.3e}')
    return errors[:10], f'{len(duo_lines)} Duo lines compared, largest relative difference {worst:.1e}'


TESTS = [('M1_2Sigma_caseb', check_caseb), ('M1_p_electron_cased', check_cased), ('M1_2Pi_2Sigma', check_compare)]

if __name__ == '__main__':
    duo, failed = Path(sys.argv[1]).resolve(), False
    for name, check in TESTS:
        with tempfile.TemporaryDirectory() as tmp:
            zpe, duo_lines = run_duo(duo, name, Path(tmp))
            strengths = levels = None
            if name in MODELS:
                levels, basis, radial = read_levels(Path(tmp), name, zpe)
                strengths = contraction(levels, basis, radial, *MODELS[name])
            errors, summary = check(duo_lines, levels, strengths)
        print(f"{'FAIL' if errors else 'PASS'} {name}: {summary}")
        for error in errors:
            print('     ', error)
        failed |= bool(errors)
    sys.exit(1 if failed else 0)
