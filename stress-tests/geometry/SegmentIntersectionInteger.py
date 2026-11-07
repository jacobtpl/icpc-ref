#!/usr/bin/env python3
"""Independent exact geometry oracle: python3 SegmentIntersectionInteger.py [driver]."""
from fractions import Fraction
import itertools
import pathlib
import random
import subprocess
import sys
import tempfile


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def intersect(points):
    a, b, c, d = points
    ab, cd = sub(b, a), sub(d, c)
    den = cross(ab, cd)
    if den:
        t = Fraction(cross(sub(c, a), cd), den)
        u = Fraction(cross(sub(c, a), ab), den)
        if not (0 <= t <= 1 and 0 <= u <= 1):
            return []
        return [(a[0] + t * ab[0], a[1] + t * ab[1])]
    direction = ab if ab != (0, 0) else cd
    if direction == (0, 0):
        return [a] if a == c else []
    if cross(sub(c, a), direction):
        return []
    left, right = max(min(a, b), min(c, d)), min(max(a, b), max(c, d))
    return [] if left > right else [left] if left == right else [left, right]


def from_hex(s):
    mantissa, exponent = s.split('p')
    sign = -1 if mantissa.startswith('-') else 1
    whole, _, tail = mantissa.lstrip('-')[2:].partition('.')
    return sign * Fraction(int(whole + tail, 16), 16**len(tail)) * Fraction(2)**int(exponent)


def run(driver):
    rng = random.Random(9183)
    n = 10**9
    cases = list(itertools.product(list(itertools.product(range(-1, 2), repeat=2)), repeat=4))
    seeds = [
        [(0, 0), (n, n), (n, 0), (0, n)],
        [(-n, -n), (n, n), (-n, n), (n, -n)],
        [(0, 0), (n, 1), (0, n), (n-1, 0)],
        [(0, 0), (701408733, 433494437), (433494437, 267914296), (267914296, 165580141)],
        [(n-1, 0), (n, 0), (n-1, 0), (n, 0)],
        [(-n, -n), (n, n), (-n+1, -n+1), (n-1, n-1)],
        [(-n, n), (n, -n), (-n, n), (-n, n)],
        [(0, 0), (1, -n), (-1, -1), (n-1, n-2)],
    ]
    for a, b, c, d in seeds:
        for p, q in [(a, b), (b, a)]:
            for r, s in [(c, d), (d, c)]:
                cases += [[p, q, r, s], [r, s, p, q]]
    for i in range(50000):
        points = [(rng.randint(-n, n), rng.randint(-n, n)) for _ in range(4)]
        mode = i % 5
        if mode == 1:
            points[1] = points[0]
        elif mode == 2:
            points[2] = points[1]
        elif mode == 3:
            points[2:] = points[:2]
        elif mode == 4:
            dx, dy = rng.randint(-10000, 10000), rng.randint(-10000, 10000)
            x, y = rng.randint(-n//2, n//2), rng.randint(-n//2, n//2)
            points = [(x+k*dx, y+k*dy) for k in rng.choices(range(-1000, 1001), k=4)]
        cases.append(points)
    for _ in range(5000):
        dx, dy = rng.randint(-100000000, 100000000), rng.randint(-100000000, 100000000)
        near = [(x+dx, y+dy) for x, y in seeds[3]]
        cases.append(near)
        cases.append([near[1], near[0], near[3], near[2]])
    data = ''.join(' '.join(str(v) for p in points for v in p)+'\n' for points in cases)
    result = subprocess.run([str(driver), '--oracle'], input=data, text=True,
                            capture_output=True, check=True)
    rows = result.stdout.splitlines()
    assert len(rows) == len(cases), (len(rows), len(cases), result.stderr)
    limits = [Fraction(1, 4000000), Fraction(1, 10**15), Fraction(1, 10**18), Fraction(1, 10**32)]
    worst = [Fraction(0)] * 4
    counts, overflows = [0, 0, 0], 0
    for points, row in zip(cases, rows):
        expected = intersect(points)
        counts[len(expected)] += 1
        blocks = row.split('|')
        assert len(blocks) == 5, row
        for k, block in enumerate(blocks[:4]):
            fields = block.split()
            assert int(fields[0]) == len(expected), (points, row, expected)
            values = list(map(from_hex, fields[1:]))
            if k == 3:
                assert len(values) == 4*len(expected), row
                values = [a+b for a, b in zip(values[::2], values[1::2])]
            assert len(values) == 2*len(expected), row
            for actual, want in zip(values, itertools.chain.from_iterable(expected)):
                if not want:
                    assert actual == 0, (k, points, row)
                error = abs(actual-want)/(abs(want) if want else 1)
                assert error <= limits[k], (k, points, row, expected, error)
                worst[k] = max(worst[k], error)
        fits = all(-(1 << 63) <= Fraction(v).numerator < (1 << 63)
                   and Fraction(v).denominator < (1 << 63)
                   for v in itertools.chain.from_iterable(expected))
        if not fits:
            assert blocks[4] == 'overflow', (points, row, expected)
            overflows += 1
        else:
            fields = blocks[4].split()
            assert int(fields[0]) == len(expected), (points, row, expected)
            assert list(map(Fraction, fields[1:])) == list(itertools.chain.from_iterable(expected)), (points, row, expected)
    print(f'Python integer intersection oracle passed ({len(cases)} cases; counts 0/1/2: {counts}; {overflows} checked overflows)')
    print('Max relative errors (float, double, long double, quad):', [f'{float(v):.3g}' for v in worst])


if __name__ == '__main__':
    if len(sys.argv) > 1:
        run(pathlib.Path(sys.argv[1]).resolve())
    else:
        with tempfile.TemporaryDirectory() as temp:
            driver = pathlib.Path(temp)/'segment'
            subprocess.run(['g++', '-std=c++17', '-O2', '-Wall', '-Wextra', '-Wconversion',
                            str(pathlib.Path(__file__).with_suffix('.cpp')), '-o', str(driver)], check=True)
            run(driver)
