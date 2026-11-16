// Driver for Template128IO.py; also self-tests the integer operators.
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main template_main
#include "../../content/contest/template.cpp"
#undef main

void check(bool ok) { if (!ok) abort(); }
template<class T> string str(T x) { ostringstream o; o << x; return o.str(); }
template<class T> bool rd(const string& s, T& x) { istringstream i(s); return bool(i >> x); }

int main() {
#ifdef __SIZEOF_INT128__
	const u128 sign = u128(1) << 127;
	const i128 hi = i128(sign - 1), lo = -hi - 1;
	check(str(lo) == "-170141183460469231731687303715884105728");
	check(str(hi) == "170141183460469231731687303715884105727");
	check(str(~u128(0)) == "340282366920938463463374607431768211455");
	check(str(u128(0)) == "0" && str(i128(0)) == "0");
	i128 x = 42; u128 u = 42;
	const char* bad[] = {"", "+", "-", "1.0", "--1", "1a", "a", "0x1"};
	rep(j,0,8) check(!rd(bad[j], x) && x == 42 && !rd(bad[j], u) && u == 42);
	check(rd("  -00012 7", x) && x == -12 && rd("+5", u) && u == 5);
	check(rd(str(lo), x) && x == lo && rd(str(hi), x) && x == hi);
	check(rd(str(~u128(0)), u) && u == ~u128(0));
	for (ll n = -20000; n <= 20000; ++n)
		check(str(i128(n)) == to_string(n) && rd(to_string(n), x) && x == n);
	unsigned long long seed = 128;
	rep(it,0,100000) {
		unsigned long long r[3];
		rep(j,0,3) r[j] = seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
		u128 b = (u128(r[0]) << 64 | r[1]) >> ((r[2] >> 33) % 128);
		check(rd(str(b), u) && u == b);
		check(rd(str(i128(b)), x) && x == i128(b));
	}
	{ ostringstream o; o << setw(5) << i128(-7) << setw(4) << u128(3); check(o.str() == "   -7   3"); }
	{ istringstream i("1 2 3"); int c = 0; while (i >> x) c++; check(c == 3); }
#endif
#if defined(__SIZEOF_INT128__) && defined(__SIZEOF_FLOAT128__)
	// "r token" -> bits of the parsed value (or "fail");
	// "w hi lo prec" -> the value with those bits, printed.
	string op;
	while (cin >> op) {
		f128 v = 0; u128 bits = 0;
		if (op == "r") {
			string t; cin >> t;
			if (!rd(t, v)) { cout << "fail\n"; continue; }
			memcpy(&bits, &v, 16);
			cout << (unsigned long long)(bits >> 64) << ' ' << (unsigned long long)bits << '\n';
		} else {
			unsigned long long h, l; int p; cin >> h >> l >> p;
			bits = u128(h) << 64 | l;
			memcpy(&v, &bits, 16);
			cout << fixed << setprecision(p) << v << '\n';
		}
	}
	cout << "f128\n";
#else
	f128 v = 0; check(rd("1.5", v) && str(v) == "1.5");
	cout << "nof128 " << numeric_limits<f128>::digits << "\n";
#endif
}
