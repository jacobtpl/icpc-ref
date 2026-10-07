#include "../utilities/template.h"

#include "../../content/numerical/MegaBigInt.h"

// Independent oracle: sign + decimal digit vector (least significant first), schoolbook.
struct Nv {
	bool neg = 0; vi d; // no leading zeros; empty = 0
	void norm() { while (!d.empty() && !d.back()) d.pop_back(); if (d.empty()) neg = 0; }
};
Nv fromStr(const string& s) {
	Nv r; int p = 0;
	if (s[0] == '-') r.neg = 1, p = 1;
	for (int i = sz(s) - 1; i >= p; i--) r.d.push_back(s[i] - '0');
	r.norm(); return r;
}
string toStr(const Nv& a) {
	if (a.d.empty()) return "0";
	string s = a.neg ? "-" : "";
	for (int i = sz(a.d); i--;) s += char('0' + a.d[i]);
	return s;
}
int cmpAbs(const Nv& a, const Nv& b) {
	if (sz(a.d) != sz(b.d)) return sz(a.d) < sz(b.d) ? -1 : 1;
	for (int i = sz(a.d); i--;) if (a.d[i] != b.d[i]) return a.d[i] < b.d[i] ? -1 : 1;
	return 0;
}
int cmp(const Nv& a, const Nv& b) {
	if (a.neg != b.neg) return a.neg ? -1 : 1;
	return a.neg ? -cmpAbs(a, b) : cmpAbs(a, b);
}
Nv addAbs(const Nv& a, const Nv& b) {
	Nv r; int c = 0;
	rep(i,0,max(sz(a.d), sz(b.d))) {
		int x = c + (i < sz(a.d) ? a.d[i] : 0) + (i < sz(b.d) ? b.d[i] : 0);
		r.d.push_back(x % 10); c = x / 10;
	}
	if (c) r.d.push_back(c);
	return r;
}
Nv subAbs(const Nv& a, const Nv& b) { // |a| >= |b|
	Nv r; int c = 0;
	rep(i,0,sz(a.d)) {
		int x = a.d[i] - c - (i < sz(b.d) ? b.d[i] : 0);
		c = x < 0; if (c) x += 10;
		r.d.push_back(x);
	}
	r.norm(); return r;
}
Nv add(const Nv& a, const Nv& b) {
	Nv r;
	if (a.neg == b.neg) r = addAbs(a, b), r.neg = a.neg;
	else if (cmpAbs(a, b) >= 0) r = subAbs(a, b), r.neg = a.neg;
	else r = subAbs(b, a), r.neg = b.neg;
	r.norm(); return r;
}
Nv neg(Nv a) { a.neg = !a.neg; a.norm(); return a; }
Nv sub(const Nv& a, const Nv& b) { return add(a, neg(b)); }
Nv mul(const Nv& a, const Nv& b) {
	Nv r; r.d.assign(sz(a.d) + sz(b.d) + 1, 0);
	rep(i,0,sz(a.d)) rep(j,0,sz(b.d)) r.d[i + j] += a.d[i] * b.d[j];
	rep(i,0,sz(r.d) - 1) r.d[i + 1] += r.d[i] / 10, r.d[i] %= 10;
	r.neg = a.neg != b.neg; r.norm(); return r;
}

mt19937_64 rng(12345);
ll rnd(ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); }

string randNum(int maxDigits, bool allowNeg = 1) {
	int len = (int)rnd(1, maxDigits), mode = (int)rnd(0, 5);
	string s;
	rep(i,0,len) {
		if (mode == 0) s += '9';
		else if (mode == 1) s += i == 0 ? '1' : '0';
		else if (mode == 2) s += rnd(0, 3) ? '0' : char('0' + rnd(0, 9));
		else if (mode == 3) s += rnd(0, 3) ? '9' : char('0' + rnd(0, 9));
		else s += char('0' + rnd(0, 9));
	}
	size_t p = s.find_first_not_of('0');
	s = p == string::npos ? "0" : s.substr(p);
	if (allowNeg && s != "0" && rnd(0, 1)) s = "-" + s;
	return s;
}
string str(const BigInt& b) { ostringstream o; o << b; return o.str(); }
string str128(__int128 x) {
	if (x == 0) return "0";
	bool n = x < 0; if (n) x = -x;
	string s; while (x) s += char('0' + (int)(x % 10)), x /= 10;
	if (n) s += '-';
	reverse(all(s)); return s;
}
#define CHECK(got, want) do { string g_ = (got), w_ = (want); if (g_ != w_) { \
	cerr << "line " << __LINE__ << ": got " << g_.substr(0, 200) << " want " << w_.substr(0, 200) << endl; assert(0); } } while (0)

ll interesting() {
	static const ll base[] = {0, 1, 2, 3, 7, 10, 999999999, 1000000000, 1000000001, 2147483647, 2147483648LL,
		4294967296LL, 999999999999999999LL, 1000000000000000000LL, 1000000000000000001LL,
		999999999000000000LL, 999999998000000001LL, LLONG_MAX};
	int m = (int)rnd(0, 3);
	ll v;
	if (m == 0) v = base[rnd(0, 17)];
	else if (m == 1) v = rnd(0, 20);
	else if (m == 2) v = rnd(0, LLONG_MAX) >> rnd(0, 62);
	else { v = base[rnd(0, 16)] + rnd(-2, 2); if (v < 0) v = 0; }
	return rnd(0, 1) ? v : -v;
}

void testSmall() {
	rep(it,0,300000) {
		ll x = interesting(), y = interesting();
		BigInt a(x), b(y);
		__int128 X = x, Y = y;
		CHECK(str(a), str128(X));
		CHECK(str(BigInt(to_string(x))), str128(X));
		CHECK(str(a + b), str128(X + Y));
		CHECK(str(a - b), str128(X - Y));
		CHECK(str(a * b), str128(X * Y));
		CHECK(str(-a), str128(-X));
		assert((a < b) == (x < y) && (a > b) == (x > y) && (a <= b) == (x <= y));
		assert((a >= b) == (x >= y) && (a == b) == (x == y) && (a != b) == (x != y));
		{ BigInt c = a; c += b; CHECK(str(c), str128(X + Y)); }
		{ BigInt c = a; c -= b; CHECK(str(c), str128(X - Y)); }
		{ BigInt c = a; c *= b; CHECK(str(c), str128(X * Y)); }
		{ BigInt c = a; c += c; CHECK(str(c), str128(2 * X)); c -= c; CHECK(str(c), "0"); }
		// mixed with built-in integers (ll must not be truncated to int)
		CHECK(str(a * y), str128(X * Y));
		{ BigInt c = a; c *= y; CHECK(str(c), str128(X * Y)); }
		CHECK(str(a + y), str128(X + Y));
		CHECK(str(x - b), str128(X - Y));
		if (y != 0) {
			// BigInt / BigInt truncates toward zero
			CHECK(str(a / b), str128(X / Y));
			{ BigInt c = a; c /= b; CHECK(str(c), str128(X / Y)); }
		}
		if (y > 0) {
			CHECK(str(a / y), str128(X / Y));
			{ BigInt c = a; c /= y; CHECK(str(c), str128(X / Y)); }
			// BigInt % BigInt (b > 0) returns the non-negative residue
			__int128 R = ((X % Y) + Y) % Y;
			CHECK(str(a % b), str128(R));
			auto qr = divmod(a, b);
			CHECK(str(qr.second), str128(R));
#ifdef BIGINT_DIVMOD_CONSISTENT
			// q*b + r == a must hold for divmod, fails for a < 0 with b not dividing a
			CHECK(str(qr.first * b + qr.second), str128(X));
#else
			if (x >= 0) CHECK(str(qr.first * b + qr.second), str128(X));
#endif
			if (y < BASE) // % long long keeps the sign of a (C++ semantics)
				assert(a % y == x % y);
		}
		if (x >= 0) {
			ll r = (ll)sqrtl((long double)x);
			while ((__int128)r * r > x) r--;
			while ((__int128)(r + 1) * (r + 1) <= x) r++;
			CHECK(str(sqrt(a)), to_string(r));
		}
		if (x && y) {
			ll g = __gcd(llabs(x), llabs(y));
			if (x > 0 && y > 0) {
				CHECK(str(gcd(a, b)), to_string(g));
				CHECK(str(lcm(a, b)), str128(X / g * Y));
			}
		}
	}
	// limits
	CHECK(str(BigInt(LLONG_MIN)), to_string(LLONG_MIN));
	CHECK(str(BigInt(LLONG_MIN) + BigInt(LLONG_MAX)), "-1");
	CHECK(str(BigInt(1) * 10000000000LL), "10000000000");
	CHECK(str(BigInt(123456789012345LL) * -4000000000LL), "-493827156049380000000000");
	CHECK(str(BigInt("100000000000000000000") / 10000000000LL), "10000000000");
	CHECK(str(BigInt("-0")), "0"); CHECK(str(BigInt("+12")), "12"); CHECK(str(BigInt("000")), "0");
	CHECK(str(BigInt("-000123")), "-123"); CHECK(str(BigInt()), "0");
	assert(BigInt("-0") == BigInt(0) && !(BigInt("-0") < BigInt(0)));
	assert(BigInt(0) * BigInt(-5) == BigInt(0) && (BigInt(5) - BigInt(5)).sign == 1);
	{ istringstream in("-123456789012345678901234567890 42"); BigInt p, q; in >> p >> q;
	  CHECK(str(p), "-123456789012345678901234567890"); CHECK(str(q), "42"); }
}

void testBig(int iters, int maxDigits) {
	rep(it,0,iters) {
		string sa = randNum(maxDigits), sb = randNum(rnd(0, 2) ? maxDigits : max(1, maxDigits / 4));
		if (rnd(0, 9) == 0) sb = sa;
		BigInt a(sa), b(sb);
		Nv A = fromStr(sa), B = fromStr(sb);
		CHECK(str(a), sa);
		CHECK(str(a + b), toStr(add(A, B)));
		CHECK(str(a - b), toStr(sub(A, B)));
		CHECK(str(a * b), toStr(mul(A, B)));
		CHECK(str(a.mul_simple(b)), toStr(mul(A, B)));
		CHECK(str(a.mul_karatsuba(b)), toStr(mul(A, B)));
		CHECK(str(a.mul_fft(b)), toStr(mul(A, B)));
		CHECK(str(a.abs()), toStr(A).substr(A.neg));
		int c = cmp(A, B);
		assert((a < b) == (c < 0) && (a == b) == (c == 0) && (a > b) == (c > 0));
		assert((a <= b) == (c <= 0) && (a >= b) == (c >= 0) && (a != b) == (c != 0));
		if (!B.d.empty()) {
			// truncated quotient: |a| = |q||b| + r', 0 <= r' < |b|, sign(q) = sign(a)sign(b)
			BigInt q = a / b;
			Nv Q = fromStr(str(q)), absA = A, absB = B, absQ = Q;
			absA.neg = absB.neg = absQ.neg = 0;
			Nv R = sub(absA, mul(absQ, absB));
			assert(!R.neg && cmpAbs(R, absB) < 0);
			assert(Q.d.empty() || Q.neg == (A.neg != B.neg));
			if (!B.neg) {
				BigInt r = a % b;
				Nv want = (A.neg && !R.d.empty()) ? sub(absB, R) : R;
				CHECK(str(r), toStr(want));
			}
			if (!A.neg && !B.neg) {
				auto qr = divmod(a, b);
				CHECK(str(qr.first), toStr(Q)); CHECK(str(qr.second), toStr(R));
			}
		}
		if (!A.neg) {
			BigInt r = sqrt(a);
			Nv Rr = fromStr(str(r)), one = fromStr("1");
			assert(!Rr.neg && cmp(mul(Rr, Rr), A) <= 0);
			Nv R1 = add(Rr, one);
			assert(cmp(mul(R1, R1), A) > 0);
		}
		int v = (int)rnd(1, rnd(0, 1) ? 20 : INT_MAX);
		{
			Nv V = fromStr(to_string(v));
			CHECK(str(a * v), toStr(mul(A, V))); CHECK(str(a * -v), toStr(mul(A, neg(V))));
			BigInt q = a / v; Nv Q = fromStr(str(q)), absA = A, absQ = Q; absA.neg = absQ.neg = 0;
			Nv R = sub(absA, mul(absQ, V));
			assert(!R.neg && cmpAbs(R, V) < 0 && (Q.d.empty() || Q.neg == A.neg));
			if (v < BASE) {
				ll m = a % (ll)v;
				CHECK(to_string(A.neg ? -m : m), toStr(R));
			}
		}
	}
}

// perfect squares and neighbours
void testSqrt() {
	rep(it,0,3000) {
		string s = randNum(80, 0);
		BigInt r(s), sq = r * r;
		CHECK(str(sqrt(sq)), s);
		CHECK(str(sqrt(sq + 1)), s == "0" ? "1" : s);
		if (s != "0") CHECK(str(sqrt(sq - 1)), str(r - 1));
	}
}

const ll P[3] = {998244353, 1000000007, 2147483629};
ll modOf(const BigInt& x, ll p) {
	ll r = 0;
	for (int i = sz(x.a); i--;) r = (ll)(((__int128)r * BASE + x.a[i]) % p);
	return x.sign < 0 ? (p - r) % p : r;
}
BigInt randLimbs(int n, int mode) {
	BigInt x; x.a.resize(n);
	for (auto& v : x.a) v = mode == 0 ? BASE - 1 : mode == 1 ? (int)rnd(0, BASE - 1) : (rnd(0, 1) ? BASE - 1 : 0);
	if (!x.a.back()) x.a.back() = 1;
	return x;
}
void checkMul(const BigInt& a, const BigInt& b, const BigInt& c) {
	assert(c.a.empty() || c.a.back() != 0);
	for (int v : c.a) assert(0 <= v && v < BASE);
	rep(k,0,3) assert(modOf(c, P[k]) == (ll)((__int128)modOf(a, P[k]) * modOf(b, P[k]) % P[k]));
	// leading limbs: product size is sz(a)+sz(b) or one less
	assert(sz(c.a) == sz(a.a) + sz(b.a) || sz(c.a) == sz(a.a) + sz(b.a) - 1);
}
void testHugeMul() {
	// operator* dispatch: simple (<= 1000111 limb products), karatsuba, fft (> 500111 limbs)
	rep(mode,0,3) {
		for (int n : {1001, 1500, 4096, 30000}) {
			BigInt a = randLimbs(n, mode), b = randLimbs(n + (int)rnd(0, 50), mode);
			if (mode == 1) b.sign = -1;
			BigInt c = a * b;
			assert(c.sign == a.sign * b.sign);
			checkMul(a, b, c);
		}
		{ // very unbalanced operands (must not be padded to equal size)
			BigInt a = randLimbs(400000, mode), b = randLimbs(3, mode);
			checkMul(a, b, a * b);
		}
		{ // FFT path through operator*
			BigInt a = randLimbs(500112, mode), b = randLimbs(1001, mode);
			checkMul(a, b, a * b);
		}
		{ // large balanced FFT
			BigInt a = randLimbs(100000, mode), b = randLimbs(100000, mode);
			checkMul(a, b, a.mul_fft(b));
		}
	}
#ifdef BIGINT_FFT_MAX
	rep(mode,0,3) {
		BigInt a = randLimbs(500112, mode), b = randLimbs(500112, mode);
		checkMul(a, b, a * b);
	}
#endif
	// division / sqrt on a few thousand limbs, verified by multiplication
	rep(mode,0,3) {
		BigInt a = randLimbs(3000, mode), b = randLimbs(1200, mode);
		auto qr = divmod(a, b);
		assert(qr.second >= 0 && qr.second < b);
		BigInt back = qr.first * b + qr.second;
		assert(back == a);
		BigInt r = sqrt(a);
		assert(r * r <= a && (r + 1) * (r + 1) > a);
	}
}

int main() {
	testSmall();
	testBig(20000, 30);
	testBig(3000, 200);
	testBig(60, 2500);
	testSqrt();
	testHugeMul();
	cout << "Tests passed!" << endl;
}
