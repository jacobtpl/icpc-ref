#include "../utilities/template.h"
#include "../../content/geometry/LineIntersectionArea.h"
#include "../../content/geometry/lineIntersection.h"
#include "../../content/number-theory/Fraction.h"

using Z = Point<ll>;
using Q = Point<Fraction>;

void check(bool ok) {
	if (!ok) abort();
}

void close(long double a, long double b) {
	check(abs(a-b) <= 8*numeric_limits<long double>::epsilon()
		* max(1.L,abs(b)));
}

int main(int argc, char**) {
	if (argc > 1) {
		Z a,b,c,d;
		cout << hexfloat;
		while (cin >> a.x >> a.y >> b.x >> b.y
				>> c.x >> c.y >> d.x >> d.y) {
			auto x=lineInterArea2<__float128>(a,b,c,d);
			long double hi=(long double)x;
			long double lo=(long double)(x-(__float128)hi);
			cout << lineInterArea2(a,b,c,d) << ' '
				<< hi << ' ' << lo << '\n';
		}
		return 0;
	}
	check(lineInterArea2(Z(0,0),Z(14,4),Z(9,0),Z(3,6))==18);
	const ll N=1000000000;
	check(lineInterArea2(Z(-N,-N),Z(N,N),Z(N,-N),Z(-N,N))==2e18L);
	check(lineInterArea2(Point<int>(-N,-N),Point<int>(N,N),
		Point<int>(N,-N),Point<int>(-N,N))==2e18L);
	using W=Point<__int128_t>;
	check(lineInterArea2<__float128>(W(-N,-N),W(N,N),
		W(N,-N),W(-N,N))==__float128(2)*N*N);
	check(lineInterArea2(Z(0,0),Z(N,1),Z(0,0),Z(1,N))==0);
	check(lineInterArea2(Z(0,0),Z(701408733,433494437),
		Z(433494437,267914296),Z(267914296,165580141))==0.5L);
	mt19937 rng(87317);
	auto coord = [&]() { return ll(rng()%61)-30; };
	auto q = [](Z p) { return Q(p.x,p.y); };
	int tested=0;
	while (tested < 20000) {
		Z a(coord(),coord()), b(coord(),coord());
		Z c(coord(),coord()), d(coord(),coord());
		if ((b-a).cross(d-c)==0) continue;
		auto intersection=lineInter(q(a),q(b),q(c),q(d));
		check(intersection.first==1);
		long double expected=(long double)q(a).cross(q(c),intersection.second);
		close(lineInterArea2(a,b,c,d),expected);
		close((long double)lineInterArea2<__float128>(a,b,c,d),expected);
		close(lineInterArea2(c,d,a,b),-expected);
		check(abs(lineInterArea2<double>(a,b,c,d)-(double)expected)
			<= 8*numeric_limits<double>::epsilon()*max(1.,abs((double)expected)));
		Z offset(coord()*1000000,coord()*1000000);
		ll scale=10000000;
		close(lineInterArea2(a*scale+offset,b*scale+offset,
			c*scale+offset,d*scale+offset),expected*scale*scale);
		++tested;
	}
	cout << "Intersection area tests passed (20000 cases)\n";
}
