/**
 * Author: Siyong
 * Date: April 6, 2024
 * Description: Modular integers for $1 \le M \le \text{INT\_MAX}$.
 *   Keep v in [0,M). Division requires a divisor coprime to M;
 *   pow requires a nonnegative integer exponent. Stream input uses ll.
 * Usage: using mint = ModInt<998244353>;
 * Status: Stress-tested
 * Time: O(1), except O(\log M) division and O(\log b) power.
 */
int const MOD = 998244353;
ll euclid(ll a, ll b, ll &x, ll &y) {
	if (!b) return x = 1, y = 0, a;
	ll d = euclid(b, a % b, y, x);
	return y -= a/b * x, d;
}
template<int M> struct ModInt {
	static_assert(M > 0);
	int v;
	constexpr explicit operator int() const {return v;}
	static constexpr int mod() {return M;}
	constexpr int val() const {return v;}
	constexpr ModInt(): v(0) {}
	constexpr ModInt(auto z) requires requires { z % M; }
		: v(int(z % M)) {
		if (v < 0) v += M;
	}
	static constexpr ModInt raw(int x) {
		assert(0 <= x && x < M);
		ModInt a; a.v = x; return a;
	}
	friend ModInt invert(ModInt a) {
		ll x, y, g = euclid(a.v, M, x, y);
		assert(g == 1); return ModInt(x);
	}
	constexpr ModInt& operator+= (ModInt const& o) {
		ll s = (ll)v + o.v;
		v = int(s >= M ? s-M : s); return *this;}
	constexpr ModInt& operator-= (ModInt const& o) {
		if ((v -= o.v) < 0) v += M;
		return *this;}
	constexpr ModInt& operator*= (ModInt const& o) {
		v = int((ll)v*o.v%M); return *this;}
	ModInt& operator/= (ModInt const& o) {
		return *this *= invert(o);}

	friend constexpr ModInt operator+ (ModInt a, ModInt const& b) {
		return a+=b;}
	friend constexpr ModInt operator- (ModInt a, ModInt const& b) {
		return a-=b;}
	friend constexpr ModInt operator* (ModInt a, ModInt const& b) {
		return a*=b;}
	friend ModInt operator/ (ModInt const& a, ModInt const& b) {
		return a*invert(b);}
	friend constexpr bool operator== (ModInt a, ModInt b) {
		return a.v == b.v;}
	friend constexpr bool operator!= (ModInt a, ModInt b) {
		return !(a == b);}

	constexpr ModInt operator+ () const {return *this;}
	constexpr ModInt operator- () const {return ModInt(-v);}
	constexpr ModInt& operator++ () {return *this += 1;}
	constexpr ModInt& operator-- () {return *this -= 1;}
	constexpr ModInt operator++ (int) {auto a = *this; ++*this; return a;}
	constexpr ModInt operator-- (int) {auto a = *this; --*this; return a;}
	friend constexpr ModInt pow(ModInt a, auto b) {
		assert(b >= 0);
		ModInt r(1);
		for(;b;b>>=1, a*=a)
			if(b&1)
				r *= a;
		return r;
	}
	friend istream& operator>> (istream& in, ModInt& a) {
		ll x; if (in >> x) a = ModInt(x);
		return in;
	}
	friend ostream& operator<< (ostream& out, ModInt a) {
		return out << a.v;
	}
};
using mint = ModInt<MOD>;
