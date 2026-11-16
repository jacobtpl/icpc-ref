#include <bits/stdc++.h>
using namespace std;

#define rep(i,a,b) for(int i = (a); i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
#define mp make_pair
#define pb push_back
using ll = long long;
using pii = pair<int,int>;
using vi = vector<int>;

template<typename T> bool ckmax(T &a, T const& b) {return b>a?a=b,1:0;}
template<typename T> bool ckmin(T &a, T const& b) {return b<a?a=b,1:0;}

// cin >> / cout << for __int128 and __float128. Pure arithmetic: no
// <quadmath.h>, -lquadmath, glibc *f128 or special flags needed.
//   i128 a; u128 b; f128 x; cin >> a >> b >> x;
//   cout << a << ' ' << fixed << setprecision(30) << x << '\n';
// Bad token => failbit (so while (cin >> a) works); ints wrap mod 2^128.
// f128 reads [+-]digits[.digits][e[+-]digits], uses the first 36
// significant digits, and always prints fixed with precision() decimals
// (default 6), about 33 significant digits correct.
// __float128 exists on x86 Linux gcc/clang and Homebrew g++ on macOS, but
// not with Apple clang; there f128 falls back to long double with the
// standard stream operators, so code still compiles and runs, just with
// less precision (64-bit on Apple Silicon): use g++-NN locally on a Mac.
#ifdef __SIZEOF_INT128__
using i128 = __int128_t;
using u128 = __uint128_t;
string u128str(u128 x) {
	string s;
	do s += char('0' + int(x % 10)); while (x /= 10);
	reverse(all(s));
	return s;
}
ostream& operator<<(ostream& o, u128 x) { return o << u128str(x); }
ostream& operator<<(ostream& o, i128 x) {
	return o << (x < 0 ? "-" + u128str(-u128(x)) : u128str(u128(x)));
}
template<class T> istream& read128(istream& in, T& x) {
	string s;
	if (!(in >> s)) return in;
	size_t i = s[0] == '-' || s[0] == '+';
	u128 v = 0;
	bool ok = i < s.size();
	for (; i < s.size(); i++) {
		ok &= !!isdigit((unsigned char)s[i]);
		v = v * 10 + u128(s[i] - '0');
	}
	if (ok) x = T(s[0] == '-' ? -v : v);
	else in.setstate(ios::failbit);
	return in;
}
istream& operator>>(istream& in, i128& x) { return read128(in, x); }
istream& operator>>(istream& in, u128& x) { return read128(in, x); }
#endif
#if defined(__SIZEOF_INT128__) && (defined(__SIZEOF_FLOAT128__) || \
	(defined(__APPLE__) && !defined(__clang__) && defined(__FLT128_MANT_DIG__)))
#define HAS_F128 1
using f128 = __float128;
istream& operator>>(istream& in, f128& x) {
	string s;
	if (!(in >> s)) return in;
	size_t i = s[0] == '-' || s[0] == '+', n = s.size();
	u128 m = 0;
	int e = 0, nd = 0, ee = 0;
	bool dot = 0, any = 0, eneg = 0;
	for (; i < n && (isdigit((unsigned char)s[i]) || (s[i] == '.' && !dot)); i++) {
		if (s[i] == '.') { dot = 1; continue; }
		any = 1;
		if (nd < 36) {
			m = m * 10 + u128(s[i] - '0'), e -= dot;
			nd += m > 0;
		} else e += !dot;
	}
	if (any && i < n && (s[i] | 32) == 'e') {
		any = 0, i++;
		if (i < n && (s[i] == '-' || s[i] == '+')) eneg = s[i++] == '-';
		for (; i < n && isdigit((unsigned char)s[i]); i++, any = 1)
			if (ee < 100000) ee = ee * 10 + (s[i] - '0');
	}
	if (!any || i < n) return in.setstate(ios::failbit), in;
	e += eneg ? -ee : ee;
	f128 r = f128(m);
	for (int a = abs(e), k; a; a -= k) { // one exact power if <= 1e48
		f128 b = 10, pw = 1;
		k = min(a, 4900);
		for (int j = k; j; j >>= 1, b *= b) if (j & 1) pw *= b;
		r = e < 0 ? r / pw : r * pw;
	}
	x = s[0] == '-' ? -r : r;
	return in;
}
ostream& operator<<(ostream& o, f128 x) {
	if (x != x) return o << "nan";
	string s = x < 0 ? "-" : "";
	if (x < 0) x = -x;
	if (x > 0 && x == x / 2) return o << s + "inf";
	int p = int(o.precision());
	string d;
	size_t n = 0, z = 0; // digits before the point, zeros right after it
	for (; x >= f128(1e37) && n < 5000; x /= f128(1e10)) p += 10, n += 10;
	// scale small x up, staying below 1 so the shifted digits follow the zeros
	for (; x > 0 && x < f128(1e-11) && p >= 10; x *= f128(1e10)) p -= 10, z += 10;
	for (; x > 0 && x < f128(.09) && p > 0; x *= 10) p--, z++;
	u128 ip = u128(x), one = u128(1) << 124;
	u128 f = u128((x - f128(ip)) * f128(one));
	d = u128str(ip), n += d.size(), d += string(z, '0');
	rep(i,0,p) f *= 10, d += char('0' + int(f >> 124)), f &= one - 1;
	if (f >= one / 2) {
		int i = sz(d) - 1;
		while (i >= 0 && d[size_t(i)] == '9') d[size_t(i--)] = '0';
		if (i < 0) d.insert(0, "1"), n++;
		else d[size_t(i)]++;
	}
	s += d.substr(0, n);
	if (n < d.size()) s += "." + d.substr(n);
	return o << s;
}
#else
using f128 = long double;
#endif

int main() { cin.tie(0);cin.sync_with_stdio(0); }
