#!/usr/bin/env python3
"""Independent arbitrary-precision oracle: python3 Fraction.py [driver]."""
import fractions
import itertools
import pathlib
import random
import subprocess
import sys
import tempfile

LOW, HIGH = -(1 << 63), (1 << 63) - 1


def checked(p, q):
    if q == 0:
        raise ZeroDivisionError
    x = fractions.Fraction(p, q)
    if not LOW <= x.numerator <= HIGH or x.denominator > HIGH:
        raise OverflowError
    return x


def expected(op, a, b, c, d):
    try:
        x, y = checked(a, b), checked(c, d)
        if op == "cmp":
            return str((x > y) - (x < y))
        if op == "ctor":
            z = x
        elif op == "neg":
            z = -x
        elif op == "abs":
            z = abs(x)
        elif op == "+":
            z = x + y
        elif op == "-":
            z = x - y
        elif op == "*":
            z = x * y
        else:
            z = x / y
        z = checked(z.numerator, z.denominator)
        return f"{z.numerator}/{z.denominator}"
    except ZeroDivisionError:
        return "domain"
    except OverflowError:
        return "overflow"


def run(driver):
    rng = random.Random(0xF12AC710)
    edges = [LOW, LOW + 1, LOW + 2, -HIGH // 2, -3, -2, -1, 0,
             1, 2, 3, HIGH // 2, HIGH - 2, HIGH - 1, HIGH]
    pairs = list(itertools.product(edges, repeat=2))
    rows = [("ctor", a, b, 0, 1) for a, b in pairs]
    for op in ("+", "-", "*", "/", "cmp"):
        rows.extend((op, a, b, c, d) for (a, b), (c, d)
                    in itertools.product(pairs, repeat=2))
    for op in ("neg", "abs"):
        rows.extend((op, a, b, 0, 1) for a, b in pairs)
    for _ in range(50000):
        values = [rng.randint(LOW, HIGH) for _ in range(4)]
        rows.append((rng.choice(["+", "-", "*", "/", "cmp"]), *values))
        # Large common factors and neighboring fractions stress cancellation.
        a, b = rng.randint(1, HIGH), rng.randint(1, HIGH)
        rows.append(("*", a, b, b, a))
        rows.append(("/", -a, b, -a, b))
        rows.append(("cmp", a, a + (a < HIGH), b, b + (b < HIGH)))
    data = "".join(" ".join(map(str, row)) + "\n" for row in rows)
    output = subprocess.run([str(driver), "--oracle"], input=data,
                            text=True, capture_output=True, check=True)
    actual = output.stdout.splitlines()
    assert len(actual) == len(rows), (len(actual), len(rows), output.stderr)
    for row, result in zip(rows, actual):
        want = expected(*row)
        assert result == want, (row, result, want)
    print(f"Python Fraction oracle passed ({len(rows)} cases)")
    run_wide(driver)


def run_wide(driver):
    rng = random.Random(0x128)
    low, high = -(1 << 127), (1 << 127) - 1
    edges = [low, low+1, low//2, LOW-1, LOW, LOW+1, -1, 0,
             1, HIGH-1, HIGH, HIGH+1, high//2, high-1, high]
    rows = list(itertools.product(edges, repeat=2))
    for _ in range(25000):
        rows.append((rng.randint(low, high), rng.randint(low, high)))
        p, q = rng.randint(LOW, HIGH), rng.randint(1, HIGH)
        scale = rng.randint(1, high // max(abs(p), q))
        rows.append((p*scale, q*scale*rng.choice([-1, 1])))
    data = ''.join(' '.join(map(str, (*divmod(p, 1 << 64), *divmod(q, 1 << 64))))
                   + '\n' for p, q in rows)
    output = subprocess.run([str(driver), '--wide'], input=data,
                            text=True, capture_output=True, check=True)
    actual = output.stdout.splitlines()
    assert len(actual) == len(rows), (len(actual), len(rows), output.stderr)
    for (p, q), result in zip(rows, actual):
        want = expected('ctor', p, q, 0, 1)
        assert result == want, (p, q, result, want)
    print(f'Python wide-ratio oracle passed ({len(rows)} cases)')


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run(pathlib.Path(sys.argv[1]).resolve())
    else:
        with tempfile.TemporaryDirectory() as temp:
            driver = pathlib.Path(temp) / "fraction"
            source = pathlib.Path(__file__).with_suffix(".cpp")
            subprocess.run(["g++", "-std=c++17", "-O2", "-Wall", "-Wextra",
                            "-Wconversion", str(source), "-o", str(driver)],
                           check=True)
            run(driver)
