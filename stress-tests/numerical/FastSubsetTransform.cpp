#include "../utilities/template.h"

#include "../../content/numerical/FastSubsetTransform.h"

mt19937_64 rng(555);
int rnd(int lo, int hi) { return lo + (int)(rng() % (unsigned)(hi - lo + 1)); }

int main() {
	rep(k,0,10) {
		vi a(1 << k), b = a, c = a, target = a;
		for(auto &x: a) x = rand() % 6 - 2;
		for(auto &x: b) x = rand() % 6 - 2;
		rep(i,0,1 << k) rep(j,0,1 << k) target[i & j] += a[i] * b[j];
		FST(a, false);
		FST(b, false);
		rep(i,0,1 << k) c[i] = a[i] * b[i];
		FST(c, true);
		assert(c == target);
	}

	// Many small cases. The header as written is the AND variant.
	rep(it,0,40000) {
		int k = rnd(0, 6), n = 1 << k, lim = rnd(0, 3) ? 30 : 1000;
		vi a(n), b(n), target(n);
		for(auto &x: a) x = rnd(-lim, lim);
		for(auto &x: b) x = rnd(-lim, lim);
		rep(i,0,n) rep(j,0,n) target[i & j] += a[i] * b[j];
		vi fa = a;
		FST(fa, false);
		FST(fa, true);
		assert(fa == a); // inverse
		assert(conv(a, b) == target);
	}
	{ // empty and single element
		vi e;
		FST(e, false); FST(e, true);
		assert(e.empty() && conv(e, e).empty());
		assert(conv({-3}, {5}) == vi{-15});
	}
	{ // large: n = 2^20, the result stays far below 2^31
		int n = 1 << 20;
		vi a(n), b(n);
		for(auto &x: a) x = rnd(-3, 3);
		for(auto &x: b) x = rnd(-3, 3);
		vi c = conv(a, b);
		// c[n-1] = a[n-1] b[n-1]; sum over all z of c[z] = (sum a)(sum b);
		// and for a fixed bit, the coefficients with that bit set are the
		// product of the corresponding half sums.
		assert(c[n-1] == a[n-1] * b[n-1]);
		ll sa = 0, sb = 0, sc = 0;
		rep(i,0,n) sa += a[i], sb += b[i], sc += c[i];
		assert(sa * sb == sc);
		rep(bit,0,20) {
			sa = sb = sc = 0;
			rep(i,0,n) if (i >> bit & 1) sa += a[i], sb += b[i], sc += c[i];
			assert(sa * sb == sc);
		}
		rep(it,0,30) { // direct evaluation of single coefficients
			int z = n - 1 - (1 << rnd(0, 19)) - (rnd(0, 1) << rnd(0, 19));
			z &= n - 1;
			// supersets x, y of z with x & y == z
			int free = (n - 1) ^ z;
			vi bits;
			rep(i,0,20) if (free >> i & 1) bits.push_back(i);
			ll want = 0;
			int m = sz(bits), tot = 1;
			rep(i,0,m) tot *= 3;
			rep(code,0,tot) {
				int x = z, y = z, t = code;
				rep(i,0,m) {
					if (t % 3 == 1) x |= 1 << bits[i];
					if (t % 3 == 2) y |= 1 << bits[i];
					t /= 3;
				}
				want += a[x] * b[y];
			}
			assert(want == c[z]);
		}
	}
	cout<<"Tests passed!"<<endl;
}
