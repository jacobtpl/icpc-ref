# Segment intersection precision

For integer input in `[-1e9,1e9]`, `segInter` separates exact intersection
classification from the scalar used to return coordinates:

```cpp
using P = Point<ll>;
P a(0,0), b(1000000000,1000000000), c(1000000000,0), d(0,1000000000);
auto usual = segInter(a,b,c,d);             // Point<long double> output
auto exact = segInter<Fraction>(a,b,c,d);   // exact if reduced coordinates fit
auto precise = segInter<__float128>(a,b,c,d);
```

| Input and output | What is guaranteed |
| --- | --- |
| Integer input, floating output | Exact vector size; approximate coordinates. The default output scalar is `long double`. |
| Integer input, `Fraction` output | Exact coordinates if their reduced numerators and denominators fit in `long long`; otherwise throws. |
| `Point<Fraction>` input | Retains `Fraction` output. Every intermediate fraction must fit, even if the answer is small. |
| Floating-point input | Retains its scalar type. Predicates and coordinates are approximate; near-collinear segments can be misclassified. |

The integer overload widens before subtraction, and computes both predicates
and complete coordinate numerators in `__int128`. It deduplicates endpoints
in exact arithmetic before converting them. Its size is always `0` (none),
`1` (one point), or `2` (a nonzero-length overlap), including degenerate point
segments. Two distinct endpoints can round to equal floating coordinates;
the vector still has size `2`. Do not deduplicate it again in floating point
if you need that classification.

The integer coordinate bound is a precondition, not a runtime check. Integer
output is rejected instead of truncating noninteger intersections. The exact
size guarantee does not extend to later predicates on rounded output points.

Changing a zero comparison to a fixed epsilon does not make the algorithm
exact. It changes the geometry being accepted and introduces a scale-dependent
tolerance. The current generic routine does not apply one implicitly.

## Why the input type matters

Consider the segments with endpoints:

```text
a = (0, 0)
b = (701408733, 433494437)
c = (433494437, 267914296)
d = (267914296, 165580141)
```

They cross at exactly `(701408733/2, 433494437/2)`. Their orientation
determinants include `1` and `-1`, obtained by subtracting products around
`1e17`. On the tested GNU C++ build, `Point<double>` loses those differences
and `segInter` returns **two** points, incorrectly reporting collinear overlap.
`Point<long double>` and `Point<Fraction>` return the correct single point.
The integer-input overload `segInter<double>` also returns the correct single
point: rounding its output cannot influence the classification.

For integer input bounded by `1e9`, a long-double implementation with at least
64 significand bits can represent these integer determinant calculations
exactly. This is a platform property, not a portable promise about the name
`long double`, and it does not apply to general noninteger coordinates.

## Representable answers can still have oversized intermediates

The diagonals of the square `(0,0), (1e9,0), (1e9,1e9), (0,1e9)` meet at
`(500000000,500000000)`. That answer fits in `Fraction`. The weighted coordinate
expression multiplies a coordinate by a determinant, producing intermediate
numerators around `1e27`. Starting with `Point<Fraction>` still throws because
each scalar operation must store its result in a Fraction.

With integer input, `segInter<Fraction>` keeps that numerator in `__int128`
until the final ratio is reduced, and returns the correct midpoint.

Normalizing the weights before multiplication helps this example, but it does
not remove bounded-fraction overflow in general. Other valid inputs need
unrepresentable final coordinates; [the fraction guide](fractions.md) includes
one such example.

For floating scalars, the same cubic products can overflow or underflow at
extreme scales even when the answer is representable. A tested square of side
`1e110` returns infinite coordinates with `double`; a side of `1e-110` produces
zero coordinates. These scales are outside the `1e9` integer-input use case.

## Choosing an approach

- For bounded exact rational constructions, use `Point<Fraction>` and account
  for intermediate growth.
- If only an intersection triangle's area is needed, use
  [`lineInterArea2`](intersection-area.md). It avoids coordinate construction
  and computes the integer numerator in `__int128` before approximation.
- For reliable segment classification on integer input with approximate
  output coordinates, use `segInter<T>` with an explicit floating output
  scalar, or the `long double` default. Keep the input coordinates integral.
- For arbitrary floating-point input, exact topology requires robust
  predicates (for example adaptive expansion arithmetic), not just a different
  interpolation expression.

The floating-input and Fraction-input paths are unchanged. `lineInter` also
continues to use its input scalar; it does not yet have the mixed-scalar mode.

## Validation

The independent Python oracle solves for exact segment parameters and uses
interval overlap for collinear cases. It checks classification and coordinates
for `float`, `double`, `long double`, `__float128` and `Fraction` output, including
checked overflow, reversals, zero-length segments, rounded-equal endpoints and
nearly collinear lines at the coordinate limit.

```sh
python3 stress-tests/geometry/SegmentIntersectionInteger.py
```
