#!/usr/bin/env python3
"""Test the 128-bit stream operators in content/contest/template.cpp.

Run: python3 stress-tests/various/Template128IO.py [compiler flags...]
Exact rational arithmetic is the oracle for __float128.
"""
from fractions import Fraction as F
from pathlib import Path
import os, random, subprocess, sys, tempfile

root = Path(__file__).resolve().parents[2]
cxx = os.environ.get('CXX', 'g++').split()
rng = random.Random(128)


def value(bits):
    sign, e, m = bits >> 127, bits >> 112 & 0x7fff, bits & (1 << 112) - 1
    v = F(m, 1) * F(2) ** -16494 if e == 0 else (m + (1 << 112)) * F(2) ** (e - 16495)
    return -v if sign else v


def ulp(v):
    v = abs(v)
    if v == 0: return F(2) ** -16494
    e = v.numerator.bit_length() - v.denominator.bit_length()
    if F(2) ** e > v: e -= 1
    return F(2) ** max(e - 112, -16494)


def dec(v, digits):
    """Decimal string of rational v with the given number of significant digits."""
    s, v = ('-' if v < 0 else ''), abs(v)
    e = int((v.numerator.bit_length() - v.denominator.bit_length()) * 0.30103) if v else 0
    v /= F(10) ** e
    while v and v >= 10: v /= 10; e += 1
    while v and v < 1: v *= 10; e -= 1
    n = str(int(v * 10 ** (digits - 1)))
    style = rng.randrange(3)
    if style == 0 or not -30 < e < 40:
        return f'{s}{n[0]}.{n[1:]}{rng.choice("eE")}{rng.choice(["", "+"]) if e >= 0 else ""}{e}'
    if e >= 0:
        n = n.ljust(e + 1, '0')
        return s + n[:e + 1] + ('.' + n[e + 1:] if n[e + 1:] or style == 1 else '')
    return s + rng.choice(['0.', '.']) + '0' * (-e - 1) + n


def main():
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'a.out'
        subprocess.run(cxx + sys.argv[1:] + [str(root / 'stress-tests/various/Template128IO.cpp'),
                                             '-o', str(exe)], check=True)
        lines, checks = [], []
        bad = ['', '.', '+', '-.', 'e5', '1e', '1e+', '1.2.3', '1x', 'nan', '--1', '1e5.0']
        for t in bad[1:]:
            lines.append(f'r {t}'); checks.append(('fail',))
        fixed = ['0', '-0.0', '1', '0.1', '-2.5', '1e0', '1E3', '.5', '5.', '+7.25e-2',
                 '0.000000000000000000000000000000000001', '123456789012345678901234567890.125',
                 '9' * 60, '0.' + '3' * 80, '1' + '0' * 50, '0.' + '0' * 50 + '1', '1e-40', '1e4000',
                 '1e-4000', '00012.5000']
        for t in fixed:
            lines.append(f'r {t}'); checks.append(('r', F(t), t))
        for _ in range(20000):
            e = rng.randrange(1, 0x7fff) if rng.random() < .1 else rng.randrange(16383 - 140, 16383 + 140)
            bits = rng.getrandbits(1) << 127 | e << 112 | rng.getrandbits(112)
            v = value(bits)
            t = dec(v, rng.choice([1, 5, 17, 34, 36, 40, 60]))
            lines.append(f'r {t}'); checks.append(('r', F(t), t))
            if abs(v) < 10 ** 60:
                p = rng.choice([0, 1, 6, 9, 15, 20, 30, 36, 45])
                lines.append(f'w {bits >> 64} {bits & (1 << 64) - 1} {p}'); checks.append(('w', v, p))
        for bits, want in [(0x7fff << 112, 'inf'), (0xffff << 112, '-inf'), (0x7fff8 << 108, 'nan')]:
            lines.append(f'w {bits >> 64} {bits & (1 << 64) - 1} 3'); checks.append((want,))
        out = subprocess.run([str(exe)], input='\n'.join(lines) + '\n', text=True,
                             capture_output=True, check=True).stdout.split('\n')
        if out[0] == 'nof128':
            print('integers passed; no __float128 on this target (f128 = long double)')
            return
        assert out[len(checks)] == 'f128', out[len(checks):]
        worst = F(0)
        fast = any(a in ('-Ofast', '-ffast-math') for a in sys.argv[1:])
        for c, got in zip(checks, out):
            if len(c) == 1 and c[0] != 'fail' and fast:
                continue  # -ffast-math assumes there are no NaNs or infinities
            if len(c) == 1:
                assert got == c[0], (c, got)
            elif c[0] == 'r':
                hi, lo = map(int, got.split())
                v = value(hi << 64 | lo)
                # only the first 36 significant digits are used
                err = abs(v - c[1])
                # one power of ten is exact up to 1e48; beyond that it drifts
                strict = c[1] == 0 or F(1, 10 ** 12) <= abs(c[1]) <= 10 ** 60
                tol = (F(3, 2) if strict else 64) * ulp(c[1]) + abs(c[1]) / 10 ** 35
                if abs(c[1]) > F(2) ** 16383 or abs(c[1]) < F(2) ** -16300: continue
                assert err <= tol, (c[2], float(err / ulp(c[1])))
                worst = max(worst, err / ulp(c[1]))
            else:
                _, v, p = c
                assert '.' in got if p else '.' not in got, got
                assert len(got.partition('.')[2]) == p, got
                err = abs(F(got) - v)
                tol = F(1, 2 * 10 ** p) + F(2) ** -124 + abs(v) / 10 ** 33
                assert err <= tol, (got, float(v), p)
                if v < 0 or got.startswith('-'): assert got.startswith('-') == (v < 0)
        print(f'Template128IO passed ({len(checks)} cases, worst read error {float(worst):.3f} ulp)')


if __name__ == '__main__':
    main()
