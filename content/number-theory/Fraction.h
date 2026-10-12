/**
 * Author: Jacob Teo, Devin
 * Date: 2026-10-12
 * License: CC0
 * Description: Exact finite fractions, stored as coprime p/q with q > 0.
 *  Uses long long storage and signed 128-bit intermediates. Throws
 *  overflow\_error if the reduced numerator or positive denominator
 *  does not fit, and domain\_error on division by zero. Do not modify
 *  p or q directly. Floating conversions are explicit; no epsilon is needed.
 *  Input accepts integer or p/q tokens; output is always p/q.
 *  fromWide(p,q) reduces a signed 128-bit ratio before checking storage.
 * Usage: using P = Point<Fraction>;
 *  Fraction x(1, 3); // x + 2 == Fraction(7, 3)
 * Time: O(log M) per arithmetic operation (128-bit Euclid).
 * Status: stress-tested
 */
#pragma once

struct Fraction {
	using I = __int128_t;
	ll p = 0, q = 1;
	Fraction(ll a = 0, ll b = 1) { *this = make(a, b); }
	static Fraction fromWide(I p, I q = 1) { return make(p, q); }
	template<class A, class B = int, enable_if_t<
		is_floating_point<A>::value ||
		is_floating_point<B>::value, int> = 0>
	Fraction(A, B = 1) = delete;

	explicit operator bool() const { return p != 0; }
	explicit operator double() const {
		return (double)p / (double)q; }
	explicit operator long double() const {
		return (long double)p / q; }
	Fraction operator+() const { return *this; }
	Fraction operator-() const { return make(-I(p), q); }
	friend Fraction operator+(Fraction a, Fraction b) {
		return make(I(a.p)*b.q + I(b.p)*a.q, I(a.q)*b.q); }
	friend Fraction operator-(Fraction a, Fraction b) {
		return make(I(a.p)*b.q - I(b.p)*a.q, I(a.q)*b.q); }
	friend Fraction operator*(Fraction a, Fraction b) {
		return make(I(a.p)*b.p, I(a.q)*b.q); }
	friend Fraction operator/(Fraction a, Fraction b) {
		return make(I(a.p)*b.q, I(a.q)*b.p); }
	Fraction& operator+=(Fraction b) { return *this = *this+b; }
	Fraction& operator-=(Fraction b) { return *this = *this-b; }
	Fraction& operator*=(Fraction b) { return *this = *this*b; }
	Fraction& operator/=(Fraction b) { return *this = *this/b; }
	friend bool operator==(Fraction a, Fraction b) {
		return a.p == b.p && a.q == b.q; }
	friend bool operator!=(Fraction a, Fraction b) {
		return !(a == b); }
	friend bool operator<(Fraction a, Fraction b) {
		return I(a.p)*b.q < I(b.p)*a.q; }
	friend bool operator>(Fraction a, Fraction b) { return b<a; }
	friend bool operator<=(Fraction a, Fraction b) {
		return !(b < a); }
	friend bool operator>=(Fraction a, Fraction b) {
		return !(a < b); }
	friend Fraction abs(Fraction a) { return a.p < 0 ? -a : a; }
	friend ostream& operator<<(ostream& os, Fraction a) {
		return os << a.p << '/' << a.q; }
	friend istream& operator>>(istream& is, Fraction& a) {
		string s;
		if (!(is >> s)) return is;
		istringstream in(s);
		ll p, q = 1; char c;
		if (!(in >> p) || (in >> c &&
			(c != '/' || !(in >> q))) || !in.eof())
			is.setstate(ios::failbit);
		else {
			try { a = Fraction(p, q); }
			catch (const domain_error&) { is.setstate(ios::failbit); }
			catch (const overflow_error&) { is.setstate(ios::failbit); }
		}
		return is;
	}
private:
	struct Raw {};
	Fraction(ll p, ll q, Raw) : p(p), q(q) {}
	static Fraction make(I p, I q) {
		if (!q) throw domain_error("Fraction: zero denominator");
		using U = __uint128_t;
		bool neg = (p < 0) != (q < 0);
		U n = p < 0 ? -U(p) : U(p), d = q < 0 ? -U(q) : U(q);
		U a = n, b = d;
		while (b) { U r = a % b; a = b; b = r; }
		n /= a; d /= a;
		if (n > U(LLONG_MAX) + neg || d > U(LLONG_MAX))
			throw overflow_error("Fraction: reduced value out of range");
		ll num = neg && n ? -ll(n-1)-1 : ll(n);
		return Fraction(num, (ll)d, Raw{});
	}
};
