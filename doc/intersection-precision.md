# Segment intersection precision

`segInter` uses the point's scalar type for both the intersection predicates
and coordinate construction. These are separate precision requirements.

| Scalar | What is guaranteed |
| --- | --- |
| `double` | Approximate predicates and coordinates; near-collinear segments can be classified incorrectly. |
| `long double` | More precision on platforms where the type is wider than `double`, but no general exactness guarantee. |
| `Fraction` | Exact predicates and coordinates when **every intermediate fraction** fits in `long long`; otherwise it throws. |
| `__int128` | Exact integer predicates for coordinates bounded by `1e9`; coordinate division still truncates noninteger answers. |

Changing a zero comparison to a fixed epsilon does not make the algorithm
exact. It changes the geometry being accepted and introduces a scale-dependent
tolerance. The current generic routine does not apply one implicitly.

## A classification failure within the 1e9 coordinate bound

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

For integer input bounded by `1e9`, a long-double implementation with at least
64 significand bits can represent these integer determinant calculations
exactly. This is a platform property, not a portable promise about the name
`long double`, and it does not apply to general noninteger coordinates.

## Representable answers can still have oversized intermediates

The diagonals of the square `(0,0), (1e9,0), (1e9,1e9), (0,1e9)` meet at
`(500000000,500000000)`. That answer fits in `Fraction`. The current weighted
coordinate expression nevertheless multiplies a coordinate by a determinant,
producing intermediate numerators around `1e27`, so `Point<Fraction>` throws.

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
  output coordinates, separate the two numeric roles: perform predicates in
  `__int128`, then convert the exact coordinate numerator and denominator to
  the chosen floating type. The current `segInter(P,P,P,P)` interface does not
  provide that mixed-scalar mode.
- For arbitrary floating-point input, exact topology requires robust
  predicates (for example adaptive expansion arithmetic), not just a different
  interpolation expression.

The generic `segInter` behavior is unchanged; this review does not claim that
one arithmetic rearrangement solves all of these cases.
