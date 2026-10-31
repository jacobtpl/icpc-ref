#!/usr/bin/env python3
"""Independent coordinate oracle: python3 LineIntersectionArea.py [driver]."""
from fractions import Fraction
import pathlib
import random
import subprocess
import sys
import tempfile


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def area(points):
    a, b, c, d = points
    ab, cd = sub(b, a), sub(d, c)
    det = cross(ab, cd)
    if not det:
        return None
    t = Fraction(cross(sub(c, a), cd), det)
    x = a[0] + t * ab[0], a[1] + t * ab[1]
    return cross(sub(c, a), sub(x, a))


def from_hex(s):
    mantissa, exponent = s.split('p')
    sign = -1 if mantissa.startswith('-') else 1
    mantissa = mantissa.lstrip('-')[2:]
    whole, _, tail = mantissa.partition('.')
    return sign * Fraction(int(whole + tail, 16), 16 ** len(tail)) * Fraction(2) ** int(exponent)


def run(driver):
    rng = random.Random(0xA2EA)
    n = 10**9
    seeds = [
        [(0, 0), (14, 4), (9, 0), (3, 6)],
        [(-n, -n), (n, n), (n, -n), (-n, n)],
        [(0, 0), (n, 1), (0, n), (n - 1, 0)],
        [(0, 0), (701408733, 433494437), (433494437, 267914296), (267914296, 165580141)],
        [(0, 0), (n, n-1), (1, 1), (n, n)],
        [(0, 0), (n, 1), (0, 0), (1, n)],
        [(0, 0), (n//2, n//2-1), (n//2, -n//2), (n-1, -2)],
    ]
    cases = []
    for a, b, c, d in seeds:
        for p, q in [(a, b), (b, a)]:
            for r, s in [(c, d), (d, c)]:
                cases += [[p, q, r, s], [r, s, p, q]]
    for _ in range(50000):
        radius = rng.choice([10, 1000, n])
        points = [(rng.randint(-radius, radius), rng.randint(-radius, radius)) for _ in range(4)]
        if area(points) is not None:
            cases.append(points)
    translated = []
    for points in cases[:2000]:
        dx, dy = rng.randint(-10**8, 10**8), rng.randint(-10**8, 10**8)
        moved = [(x + dx, y + dy) for x, y in points]
        if all(abs(x) <= n and abs(y) <= n for x, y in moved):
            translated.append(moved)
    cases += translated
    data = ''.join(' '.join(str(v) for p in points for v in p) + '\n' for points in cases)
    result = subprocess.run([str(driver), '--oracle'], input=data, text=True,
                            capture_output=True, check=True)
    rows = result.stdout.splitlines()
    assert len(rows) == len(cases), (len(rows), len(cases), result.stderr)
    worst_ld = worst_quad = Fraction(0)
    for points, row in zip(cases, rows):
        expected = area(points)
        ld, hi, lo = map(from_hex, row.split())
        error_ld = abs(ld - expected) / max(1, abs(expected))
        error_quad = abs(hi + lo - expected) / max(1, abs(expected))
        assert error_ld < Fraction(1, 10**18), (points, row, expected)
        assert error_quad < Fraction(1, 10**32), (points, row, expected)
        worst_ld, worst_quad = max(worst_ld, error_ld), max(worst_quad, error_quad)
    print(f'Python intersection-coordinate oracle passed ({len(cases)} cases)')
    print(f'Max normalized error: long double={float(worst_ld):.3g}, __float128={float(worst_quad):.3g}')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        run(pathlib.Path(sys.argv[1]).resolve())
    else:
        with tempfile.TemporaryDirectory() as temp:
            driver = pathlib.Path(temp) / 'area'
            subprocess.run(['g++', '-std=c++17', '-O2', '-Wall', '-Wextra', '-Wconversion',
                            str(pathlib.Path(__file__).with_suffix('.cpp')), '-o', str(driver)], check=True)
            run(driver)
