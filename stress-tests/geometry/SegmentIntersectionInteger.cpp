#include "../utilities/template.h"
#include "../../content/geometry/SegmentIntersection.h"

using Z = Point<ll>;
using Q = Point<Fraction>;

void check(bool ok) { if (!ok) abort(); }

template<class T> void emit(Z a, Z b, Z c, Z d) {
	try {
		auto out=segInter<T>(a,b,c,d);
		cout << out.size();
		for (auto p : out) for (auto v : {p.x,p.y}) {
			if constexpr (is_same<T,__float128>::value) {
				long double hi=(long double)v;
				cout << ' ' << hi << ' ' << (long double)(v-(__float128)hi);
			} else cout << ' ' << v;
		}
	} catch (const overflow_error&) { cout << "overflow"; }
}

int main(int argc, char**) {
	if (argc > 1) {
		Z a,b,c,d;
		cout << hexfloat;
		while (cin >> a.x >> a.y >> b.x >> b.y
				>> c.x >> c.y >> d.x >> d.y) {
			emit<float>(a,b,c,d); cout << '|';
			emit<double>(a,b,c,d); cout << '|';
			emit<long double>(a,b,c,d); cout << '|';
			emit<__float128>(a,b,c,d); cout << '|';
			emit<Fraction>(a,b,c,d); cout << '\n';
		}
		return 0;
	}
	static_assert(is_same<decltype(segInter(Z(),Z(),Z(),Z())),
		vector<Point<long double>>>::value);
	const ll N=1000000000;
	check((segInter<Fraction>(Z(0,0),Z(N,N),Z(N,0),Z(0,N)) ==
		vector<Q>{Q(N/2,N/2)}));
	check((segInter<Fraction>(Z(-N,-N),Z(N,N),Z(-N,N),Z(N,-N)) ==
		vector<Q>{Q(0,0)}));
	check((segInter<Fraction>(Z(0,0),Z(1,-N),Z(-1,-1),Z(N-1,N-2)) ==
		vector<Q>{Q(Fraction(1,N*N+N-1),Fraction(-N,N*N+N-1))}));
	check((segInter<Fraction>(Z(0,0),Z(701408733,433494437),
		Z(433494437,267914296),Z(267914296,165580141)) ==
		vector<Q>{Q(Fraction(701408733,2),Fraction(433494437,2))}));
	check(segInter<double>(Z(0,0),Z(701408733,433494437),
		Z(433494437,267914296),Z(267914296,165580141)).size()==1);
	auto overlap=segInter<float>(Z(N-1,0),Z(N,0),Z(N-1,0),Z(N,0));
	check(overlap.size()==2 && overlap[0]==overlap[1]);
	check((segInter<Fraction>(Z(2,2),Z(2,2),Z(0,0),Z(4,4)) ==
		vector<Q>{Q(2,2)}));
	check(segInter(Z(2,1),Z(2,1),Z(0,0),Z(4,4)).empty());
	bool overflow=false;
	try { segInter<Fraction>(Z(0,0),Z(N,1),Z(0,N),Z(N-1,0)); }
	catch (const overflow_error&) { overflow=true; }
	check(overflow);
	check(segInter<__float128>(Z(0,0),Z(N,1),Z(0,N),Z(N-1,0)).size()==1);
	check(segInter(Point<int>(-N,-N),Point<int>(N,N),
		Point<int>(N,-N),Point<int>(-N,N)).size()==1);
	using W=Point<__int128_t>;
	check((segInter<Fraction>(W(0,0),W(N,N),W(N,0),W(0,N)) ==
		vector<Q>{Q(N/2,N/2)}));
	mt19937 rng(9183);
	auto coord=[&]() { return ll(rng()%61)-30; };
	auto q=[](Z p) { return Q(p.x,p.y); };
	for (int i=0; i<20000; ++i) {
		Z a(coord(),coord()), b(coord(),coord()), c(coord(),coord()), d(coord(),coord());
		auto exact=segInter(q(a),q(b),q(c),q(d));
		check(segInter<Fraction>(a,b,c,d)==exact);
		check(segInter<Fraction>(b,a,d,c)==exact);
		check(segInter<Fraction>(c,d,a,b)==exact);
		auto approx=segInter(a,b,c,d);
		check(approx.size()==exact.size());
		for (size_t j=0; j<exact.size(); ++j) {
			long double x=(long double)exact[j].x, y=(long double)exact[j].y;
			check(abs(approx[j].x-x)<=8*numeric_limits<long double>::epsilon()*max(1.L,abs(x)));
			check(abs(approx[j].y-y)<=8*numeric_limits<long double>::epsilon()*max(1.L,abs(y)));
		}
	}
	cout << "Integer segment intersection tests passed (20000 cases)\n";
}
