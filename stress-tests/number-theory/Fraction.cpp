#include "../utilities/template.h"
#include "../../content/number-theory/Fraction.h"

using I = __int128_t;

void check(Fraction x, I p, I q) {
	assert(x.q > 0);
	assert(I(x.p)*q == p*x.q);
	I a = x.p < 0 ? -I(x.p) : x.p, b = x.q;
	while (b) { I r = a % b; a = b; b = r; }
	assert(a == 1);
}

template<class E, class F> void throws(F f) {
	bool caught = false;
	try { f(); } catch (const E&) { caught = true; }
	assert(caught);
}

int main(int argc, char** argv) {
	if (argc > 1 && string(argv[1]) == "--wide") {
		ll ph, qh;
		unsigned long long pl, ql;
		const I base = I(1)<<64;
		while (cin >> ph >> pl >> qh >> ql) {
			try {
				cout << Fraction::fromWide(I(ph)*base+I(pl), I(qh)*base+I(ql));
			} catch (const domain_error&) { cout << "domain"; }
			catch (const overflow_error&) { cout << "overflow"; }
			cout << '\n';
		}
		return 0;
	}
	if (argc > 1) {
		string op; ll a, b, c, d;
		while (cin >> op >> a >> b >> c >> d) {
			try {
				Fraction x(a,b), y(c,d);
				if (op == "ctor") cout << x;
				if (op == "+") cout << x+y;
				if (op == "-") cout << x-y;
				if (op == "*") cout << x*y;
				if (op == "/") cout << x/y;
				if (op == "neg") cout << -x;
				if (op == "abs") cout << abs(x);
				if (op == "cmp") cout << ((x>y)-(x<y));
			} catch (const domain_error&) { cout << "domain"; }
			catch (const overflow_error&) { cout << "overflow"; }
			cout << '\n';
		}
		return 0;
	}
	static_assert(!is_constructible<Fraction, double>::value);
	static_assert(!is_constructible<Fraction, int, double>::value);
	static_assert(!is_convertible<Fraction, double>::value);
	static_assert(is_convertible<ll, Fraction>::value);
	check(Fraction(), 0, 1);
	check(Fraction(LLONG_MIN, LLONG_MIN), 1, 1);
	check(Fraction(0, LLONG_MIN), 0, 1);
	check(Fraction(2, LLONG_MIN), -1, I(1)<<62);
	check(Fraction(LLONG_MIN, 2), I(LLONG_MIN)/2, 1);
	check(Fraction(LLONG_MAX, 2)*Fraction(2, LLONG_MAX), 1, 1);
	check(Fraction(LLONG_MIN)-Fraction(LLONG_MIN), 0, 1);
	check(Fraction(LLONG_MIN)/Fraction(LLONG_MIN), 1, 1);
	check(Fraction(LLONG_MIN, LLONG_MAX)+Fraction(1), -1, LLONG_MAX);
	const I lo = -(I(1)<<126)-(I(1)<<126), hi = -(lo+1);
	assert(Fraction::fromWide(lo,lo) == 1);
	assert(Fraction::fromWide(hi,hi) == 1);
	assert(Fraction::fromWide(0,lo) == 0);
	assert(Fraction::fromWide(lo,lo/2) == 2);
	assert(Fraction::fromWide(lo,I(1)<<64) == Fraction(LLONG_MIN));
	assert(Fraction::fromWide(42) == 42);
	throws<overflow_error>([&] { Fraction::fromWide(lo,-1); });
	throws<overflow_error>([&] { Fraction::fromWide(1,lo); });
	throws<overflow_error>([&] { Fraction::fromWide(hi,lo); });
	throws<domain_error>([&] { Fraction::fromWide(lo,0); });
	throws<domain_error>([] { Fraction(0,0); });
	throws<domain_error>([] { Fraction(1)/0; });
	throws<overflow_error>([] { Fraction(1, LLONG_MIN); });
	throws<overflow_error>([] { Fraction(LLONG_MIN, -1); });
	throws<overflow_error>([] { -Fraction(LLONG_MIN); });
	throws<overflow_error>([] { abs(Fraction(LLONG_MIN)); });
	throws<overflow_error>([] { Fraction(LLONG_MAX)+1; });
	throws<overflow_error>([] { Fraction(LLONG_MIN)-1; });
	throws<overflow_error>([] { Fraction(1, LLONG_MAX)/2; });
	Fraction keep(LLONG_MAX);
	throws<overflow_error>([&] { keep += 1; });
	assert(keep == LLONG_MAX);
	throws<domain_error>([&] { keep /= 0; });
	assert(keep == LLONG_MAX);
	assert((bool)Fraction(-1) && !(bool)Fraction());
	assert((double)Fraction(1,2) == 0.5);
	assert((long double)Fraction(-1,4) == -0.25L);
	for (string s : {"", "1/0", "1/", "1/2/3", "1.2", "x",
		"9223372036854775808", "1/-9223372036854775808"}) {
		Fraction x(7,8); istringstream in(s); in >> x;
		assert(in.fail() && x == Fraction(7,8));
	}
	istringstream in("0 -3 2/-4 -9223372036854775808/2");
	for (Fraction want : {Fraction(), Fraction(-3), Fraction(-1,2),
		Fraction(LLONG_MIN,2)}) {
		Fraction x; in >> x; assert(!in.fail() && x == want);
	}
	mt19937_64 rng(0xF12AC710);
	rep(it,0,200000) {
		ll a = ll(rng()%20001)-10000, b = ll(rng()%200)+1;
		ll c = ll(rng()%20001)-10000, d = ll(rng()%200)+1;
		if (rng()&1) b = -b;
		if (rng()&1) d = -d;
		Fraction x(a,b), y(c,d);
		check(x, a, b); check(y, c, d);
		check(x+y, I(a)*d+I(c)*b, I(b)*d);
		check(x-y, I(a)*d-I(c)*b, I(b)*d);
		check(x*y, I(a)*c, I(b)*d);
		if (c) check(x/y, I(a)*d, I(b)*c);
		I lhs = I(a)*d, rhs = I(c)*b;
		int cmp = ((lhs>rhs)-(lhs<rhs)) * ((b>0)==(d>0) ? 1 : -1);
		assert((x<y) == (cmp<0) && (x<=y) == (cmp<=0));
		assert((x>y) == (cmp>0) && (x>=y) == (cmp>=0));
		assert((x==y) == (cmp==0) && (x!=y) == (cmp!=0));
		assert((x+y)-y == x && x+y == y+x && x*y == y*x);
		assert(2+x == x+2 && 2*x == x*2 && (2-x)+x == 2);
		if (x != 0) assert((2/x)*x == 2);
		Fraction z=x; z+=y; assert(z==x+y);
		z=x; z-=y; assert(z==x-y);
		z=x; z*=y; assert(z==x*y);
		if (y != 0) { z=x; z/=y; assert(z==x/y); }
		z=x; z+=z; assert(z==2*x);
		z=x; z-=z; assert(z==0);
		z=x; z*=z; assert(z==x*x);
		if (x != 0) { z=x; z/=z; assert(z==1); }
		stringstream ss; ss << x; ss >> z; assert(z==x);
	}
	cout << "Fraction tests passed (200000 random pairs)\n";
}
