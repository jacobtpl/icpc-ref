/**
 * Author: Simon Lindholm
 * Date: 2016-07-24
 * License: CC0
 * Source: folklore
 * Description: Finds \{x, y, d\} s.t. ax + by = d = gcd(a, b)
 * (d = -gcd(a, b) is possible if a or b is negative).
 * Status: Tested
 */

static BigInteger[] euclid(BigInteger a, BigInteger b) {
	BigInteger x = BigInteger.ONE, yy = x;
	BigInteger y = BigInteger.ZERO, xx = y;
	while (b.signum() != 0) {
		BigInteger q = a.divide(b), t = b;
		b = a.remainder(b); a = t;
		t = xx; xx = x.subtract(q.multiply(xx)); x = t;
		t = yy; yy = y.subtract(q.multiply(yy)); y = t;
	}
	return new BigInteger[]{x, y, a};
}
