# Exact fractions

Paste `content/number-theory/Fraction.h` after the contest template.
It uses `long long` numerator/denominator storage and GNU `__int128`
intermediates. It has no external dependencies. C++17 is needed for the
numeric-type dispatch in the algebra and half-plane routines.

```cpp
Fraction a(1, 3), b(2, 5);
assert(a + b == Fraction(11, 15));
assert(2 * a == Fraction(2, 3));
cout << a;               // 1/3
cin >> b;                // one integer or p/q token, e.g. -6/8
double approximate = double(a);
```

`p` and `q` are coprime, `q > 0`, and zero is `0/1`. Do not mutate these
fields. Arithmetic, all six comparisons, compound assignments, unary signs,
`abs`, explicit `bool`, and explicit `double`/`long double` conversions are
provided. Integer literals work on either side of arithmetic and comparison.
Ordinary constructor arguments must fit in `long long`; floating-point arguments
are rejected instead of silently truncating. `Fraction::fromWide(p,q)` accepts
signed `__int128` arguments, reduces them before checking the storage limits,
and throws if the reduced result still does not fit. To represent a decimal exactly,
write `Fraction(1, 10)`, not `0.1`.

## Range and cost

Every result is reduced **before** narrowing back to `long long`. Products,
cross products used for comparison, and sums of two numerator products fit in
signed `__int128` for valid stored fractions. This includes `LLONG_MIN`:
normalization never negates it in `long long` or calls `std::gcd` on it.
Normalization uses unsigned magnitudes, so `fromWide` also handles the full
signed 128-bit range, including its minimum value.

For example, `Fraction(LLONG_MAX, 2) * Fraction(2, LLONG_MAX)` is exactly one,
and `Fraction(LLONG_MIN, LLONG_MIN)` is also one. An unrepresentable reduced
result throws `overflow_error`; a zero denominator or division by zero throws
`domain_error`. Failed compound assignment leaves its left operand unchanged.
Invalid input sets the stream's fail bit and leaves the destination unchanged.

This is a bounded exact type, not a big integer. Every intermediate **fraction**
in an algorithm must fit, even if the algorithm's final answer is small.
Repeated geometric constructions can grow denominators quickly. In particular,
`-Fraction(LLONG_MIN)` and `abs(Fraction(LLONG_MIN))` cannot be represented.
Keep overflow checks enabled; compiling with `NDEBUG` does not remove them.
Arithmetic runs Euclid on at most 128-bit integers, while comparisons are
constant-time cross products. There is no per-testcase initialization or pool.

Even integer coordinates bounded by `1e9` can produce intersection coordinates
whose reduced numerators do not fit. For example, the segments
`(0,0)--(1000000000,1)` and `(0,1000000000)--(999999999,0)` intersect at
an x-coordinate of `999999999000000000000000000/1000000000999999999`.
`segInter` with `Point<Fraction>` correctly throws on this input.

For this coordinate bound, ratios describing the position along the original
segments still fit: their numerators and denominators are determinants of
integer differences, bounded in magnitude by `8e18`. Compute determinants
before constructing the ratio, using `Point<__int128_t>` where needed.
An area-only calculation can often use these ratios without constructing the
intersection coordinates. Accumulating many rational areas still needs a
separate numeric choice: denominators can grow beyond `long long`, and
subtracting large areas can require more floating-point precision.

For triangle areas involving an intersection of integer-coordinate lines,
[`lineInterArea2`](intersection-area.md) computes exact wide integer
intermediates followed by one floating-point division, without constructing
the intersection or a `Fraction`. Its derivation and diagram explain the
shared-base area ratio.

The [segment precision review](intersection-precision.md) gives concrete
examples of floating-point misclassification and avoidable intermediate
growth, and explains how the integer-input overload avoids these failures.

## Geometry

For integer input with coordinates in `[-1e9,1e9]`, keep the input points
integral and select the result scalar when calling `segInter`:

```cpp
using Z = Point<ll>;
Z a(0,0), b(1000000000,1000000000), c(1000000000,0), d(0,1000000000);
auto exact = segInter<Fraction>(a,b,c,d); // (500000000,500000000)
auto approx = segInter<__float128>(a,b,c,d);
auto usual = segInter(a,b,c,d);          // vector<Point<long double>>
```

This overload computes predicates and the complete coordinate numerators in
`__int128`. It converts only the finished answers, so a representable Fraction
coordinate is not rejected because a temporary numerator exceeded `long long`.
The size is exact even with floating output: `0` for none, `1` for one point,
and `2` for a nonzero-length overlap. Deduplication happens before conversion;
two distinct overlap endpoints may round to equal floating points while the
vector still correctly has size `2`. Integer output is not supported.

The bound is a precondition, not a runtime range check. For noninteger inputs,
the existing `segInter(Point<Fraction>,...)` and floating-input calls retain
their scalar type and its arithmetic limits. The integer overload does not
change `lineInter`, which still uses its input scalar for all arithmetic.

```cpp
using Q = Point<Fraction>;
auto [kind, p] = lineInter(Q(0,0), Q(3,3), Q(0,2), Q(3,0));
assert(kind == 1 && p == Q(Fraction(6,5), Fraction(6,5)));
vector<Q> triangle = {Q(0,0), Q(2,0), Q(0,1)};
auto h = convexHull(triangle);
Fraction area = polygonArea2(h) / 2;
```

The following accept fraction points directly, with no epsilon:

- `Point` / `Point3D`: arithmetic, ordering, equality, dot/cross products,
  squared distances, and 2D perpendicular vectors.
- `sideOf` (the overload without epsilon), `onSegment`, `lineInter`,
  `segInter`, `lineProj` (projection and reflection).
- `convexHull`, `hullDiameter`, `inHull`, `inPolygon`, `lineHull`, `extrVertex`.
- `polygonArea2`, `polygonCenter`, `polygonCut`, `polyUnion`,
  `linearTransformation`, `ccCenter`, and `signedPolyVolume`.
- `halfPlaneIsect` with `vector<RayT<Fraction>>`; the existing `Ray` alias is
  still `RayT<double>`. Fraction instantiations use exactly zero tolerance.
  Existing preconditions still apply: bounded, nonzero-area intersection;
  valid nonzero ray directions. `add_bounds` clips to the existing
  `[0, 10^9]` square and may exceed fraction storage limits in constructions.

Notebook headers retain their existing example `P`/`T` aliases for existing
callers. When combining snippets, choose one alias (e.g. `using P =
Point<Fraction>`) rather than pasting conflicting `Point<ll>` and
`Point<double>` aliases. The templates also accept `Point<Fraction>` under
another name without changing the example aliases. Existing geometry
preconditions, such as a nonempty hull or a nondegenerate triangle, remain.

Lengths, normalized vectors, arbitrary rotations, and general circle
intersection coordinates are not necessarily rational. `dist`, `angle`,
`phi`, and `theta` return approximate doubles. For `unit`, `normal`, `rotate`,
circle routines, and other floating-point-only snippets, explicitly construct
`Point<double>(double(p.x), double(p.y))` (and similarly in 3D). Those results
do not have exact predicates. Use `dist2` for exact distance comparisons.

## Algebra and other snippets

`Matrix<Fraction, N>`, `det`, `gauss`, `matInv`, `solveLinear`, `tridiagonal`,
and `interpolate` accept fraction data directly. Elimination uses exact zero
tests for fractions and preserves the existing tolerances for floating types.
For example, a pivot of `Fraction(1, 1000000000000LL)` is nonzero.

Other notebook snippets with a manually selected scalar alias can use
`Fraction` when they need only field arithmetic and comparisons. Replace
floating constants with exact fractions and epsilon tests with zero tests.
For example, a scalar `max(x, 0.0)` becomes `max(x, Fraction(0))`.

This does not turn every snippet into a rational algorithm: `Simplex` still
uses floating infinities as result sentinels, and needs a separate status
representation to use finite fractions. Integer/modular algorithms (NTT,
modular FPS, bit operations, integer grid optimizations), transcendental
functions, and numerical integration/root finding retain their original
domains. There is no implicit conversion that silently approximates a
fraction to make these compile.

## Reproduce the tests

From the repository root:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -Wconversion stress-tests/number-theory/Fraction.cpp -o /tmp/fraction
/tmp/fraction
python3 stress-tests/number-theory/Fraction.py /tmp/fraction
g++ -std=c++17 -O2 -Wall -Wextra -Wconversion stress-tests/geometry/Fraction.cpp -o /tmp/fraction-geometry
/tmp/fraction-geometry
g++ -std=c++17 -O2 -Wall -Wextra -Wconversion stress-tests/numerical/Fraction.cpp -o /tmp/fraction-algebra
/tmp/fraction-algebra
```

The Python test uses the standard library's arbitrary-precision
`fractions.Fraction` as an independent oracle, including normalization,
`LLONG_MIN`/`LLONG_MAX`, cancellation, ordering, and overflow classification.
The C++ suites exercise arithmetic identities, failure behavior, random exact
geometry and linear algebra. They are also included in `make test`.
