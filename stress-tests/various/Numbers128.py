#!/usr/bin/env python3
"""Test integer I/O and the quadmath program printed in the reference.

Run: python3 stress-tests/various/Numbers128.py [--sanitize]
Requires g++, libquadmath, and Python's standard library.
"""
import argparse
from decimal import Decimal, localcontext
from pathlib import Path
import random
import re
import subprocess
import tempfile


def integers(driver):
    low, high = -(1 << 127), (1 << 127) - 1
    rng = random.Random(128)
    tokens = ['', '+', '-', ' 1', '1 ', '1.0', '++1', '--1', '+-1',
              '0x12', '12a', '-0', '+000', '9' * 1000, '0' * 1000]
    for n in [low-1, low, low+1, -(1 << 63), -1, 0, 1,
              (1 << 63)-1, high-1, high, high+1]:
        tokens += [str(n), ('-' if n < 0 else '+') + '000' + str(abs(n))]
    tokens += [str(rng.randint(low*2, high*2)) for _ in range(20000)]
    result = subprocess.run([str(driver), '--oracle'],
                            input='\n'.join(tokens)+'\n', text=True,
                            capture_output=True, check=True)
    output = result.stdout.splitlines()
    assert len(output) == len(tokens)
    for token, actual in zip(tokens, output):
        value = int(token) if re.fullmatch(r'[+-]?[0-9]+', token) else None
        expected = str(value) if value is not None and low <= value <= high else 'invalid'
        assert actual == expected, (token, actual, expected)
    print(f'Python integer oracle passed ({len(tokens)} cases)')


def quadmath(driver):
    def convert(token):
        return subprocess.run([str(driver)], input=token+'\n', text=True,
                              capture_output=True, check=True).stdout.strip()

    rng = random.Random(113)
    tokens = ['0', '-0', '1', '-1', '0.1', '1e4932', '-1e4932',
              '1e-4930', '-1e-4930',
              '1.0000000000000000000000000000000002']
    tokens += [f'{rng.randint(-(1 << 127), (1 << 127)-1)}e{rng.randint(-4900, 4800)}'
               for _ in range(100)]
    with localcontext() as ctx:
        ctx.prec = 100
        tolerance = Decimal(2) ** -112
        for token in tokens:
            output = convert(token)
            assert convert(output) == output, (token, output)
            expected, actual = Decimal(token), Decimal(output)
            if expected:
                assert abs((actual-expected)/expected) < tolerance, (token, output)
            else:
                assert actual == 0
                assert actual.is_signed() == expected.is_signed()
        assert Decimal(convert(tokens[9])) > 1
        tiny = convert(str(Decimal(2) ** -16494))
        assert Decimal(tiny) > 0 and convert(tiny) == tiny
    print(f'Quadmath reference example passed ({len(tokens)+1} cases and round trips)')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sanitize', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    reference = (root/'content/various/Numbers128.tex').read_text()
    match = re.search(r'\\begin\{lstlisting\}\[language=C\+\+,label=quadmath-example\]\n(.*?)\\end\{lstlisting\}',
                      reference, re.S)
    assert match, 'Quadmath example missing from reference'
    flags = ['g++', '-std=c++17', '-O2', '-Wall', '-Wextra', '-Wconversion', '-Werror']
    if args.sanitize:
        flags += ['-O1', '-g', '-fsanitize=address,undefined',
                  '-fno-sanitize-recover=all', '-fno-omit-frame-pointer']
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        int_driver, quad_driver = temp/'integer', temp/'quad'
        source = temp/'quad.cpp'
        source.write_text(match[1])
        subprocess.run(flags + [str(root/'stress-tests/various/Int128IO.cpp'),
                                '-o', str(int_driver)], check=True)
        subprocess.run(flags + [str(source), '-o', str(quad_driver), '-lquadmath'], check=True)
        subprocess.run([str(int_driver)], check=True)
        integers(int_driver)
        quadmath(quad_driver)


if __name__ == '__main__':
    main()
