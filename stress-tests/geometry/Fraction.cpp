#include "../utilities/template.h"
#include "../../content/number-theory/Fraction.h"
#include "../../content/geometry/Point.h"
#include "../../content/geometry/Point3D.h"
#include "../../content/geometry/sideOf.h"
#include "../../content/geometry/PolygonArea.h"
#include "../../content/geometry/lineIntersection.h"
#include "../../content/geometry/SegmentIntersection.h"
#include "../../content/geometry/LineProjRefl.h"
#include "../../content/geometry/LineHullIntersection.h"
#include "../../content/geometry/PolyhedronVolume.h"

namespace hull {
#include "../../content/geometry/ConvexHull.h"
#include "../../content/geometry/HullDiameter.h"
#include "../../content/geometry/PointInsideHull.h"
}
namespace polygon {
#include "../../content/geometry/SegmentDistance.h"
#include "../../content/geometry/InsidePolygon.h"
#include "../../content/geometry/PolygonCut.h"
#include "../../content/geometry/PolygonCenter.h"
#include "../../content/geometry/PolygonUnion.h"
#include "../../content/geometry/linearTransformation.h"
#include "../../content/geometry/circumcircle.h"
}
namespace halfplane {
#define pb push_back
#include "../../content/geometry/HalfPlaneIntersection.h"
#undef pb
#undef reduce
}

using Q = Point<Fraction>;
using Z = Point<ll>;

vector<Q> segmentOracle(Z a, Z b, Z c, Z d) {
	ll den = (b-a).cross(d-c);
	auto q = [](Z p) { return Q(p.x, p.y); };
	if (den) {
		Fraction t((c-a).cross(d-c), den);
		Fraction u((c-a).cross(b-a), den);
		if (t<0 || t>1 || u<0 || u>1) return {};
		return {q(a)+(q(b)-q(a))*t};
	}
	vector<Q> out;
	for (Z p : {a,b,c,d}) {
		auto contains = [](Z x, Z y, Z z) {
			return (y-x).cross(z-x)==0 &&
				min(x.x,y.x)<=z.x && z.x<=max(x.x,y.x) &&
				min(x.y,y.y)<=z.y && z.y<=max(x.y,y.y);
		};
		if (contains(a,b,p) && contains(c,d,p)) out.push_back(q(p));
	}
	sort(all(out)); out.erase(unique(all(out)), out.end());
	return out;
}

Fraction area(vector<Q> p) {
	return p.empty() ? Fraction(0) : abs(polygonArea2(p))/2;
}

vector<Q> rect(Fraction x, Fraction y, Fraction w, Fraction h) {
	return {Q(x,y), Q(x+w,y), Q(x+w,y+h), Q(x,y+h)};
}

int main() {
	Q zero, a(1, Fraction(1,1000000000));
	Q b(1, Fraction(1,1000000001));
	assert(sideOf(zero,a,b)==-1);
	assert(!onSegment(zero,a,b));
	Q big(1000000000,1000000000);
	assert(sideOf(big,big+Q(1,0),big+Q(0,1))==1);
	assert(lineInter(Q(0,0),Q(3,3),Q(0,2),Q(3,0)).second ==
		Q(Fraction(6,5),Fraction(6,5)));
	assert(lineInter(Q(0,0),Q(1,1),Q(2,2),Q(3,3)).first==-1);
	assert(lineInter(Q(0,0),Q(1,1),Q(0,1),Q(1,2)).first==0);
	assert((segInter(Q(0,0),Q(3,3),Q(0,2),Q(3,0))==
		vector<Q>{Q(Fraction(6,5),Fraction(6,5))}));
	assert((segInter(Q(3,0),Q(0,0),Q(1,0),Q(4,0))==
		vector<Q>{Q(1,0),Q(3,0)}));
	assert((segInter(Q(0,0),Q(1,1),Q(1,1),Q(2,0))==vector<Q>{Q(1,1)}));
	assert((segInter(Q(1,1),Q(1,1),Q(0,0),Q(2,2))==vector<Q>{Q(1,1)}));
	assert(segInter(Q(0,0),Q(1,0),Q(2,0),Q(3,0)).empty());
	const ll N=1000000000;
	assert((segInter(Q(0,0),Q(1,Fraction(1,N)),
		Q(0,Fraction(1,N)),Q(1,Fraction(1,N+1)))==
		vector<Q>{Q(Fraction(N+1,N+2),Fraction(N+1,N*(N+2)))}));
	auto small=rect(0,0,Fraction(1,1000000),Fraction(1,1000000));
	vector<halfplane::RayT<Fraction>> tinyRays;
	rep(i,0,4) tinyRays.push_back({small[i],small[(i+1)%4]-small[i]});
	tinyRays.push_back(tinyRays[0]);
	assert(area(halfplane::halfPlaneIsect(tinyRays))==Fraction(1,1000000000000LL));
	mt19937 rng(127);
	auto coord = [&]() { return int(rng()%17)-8; };
	rep(it,0,20000) {
		Z za(coord(),coord()), zb(coord(),coord());
		Z zc(coord(),coord()), zd(coord(),coord());
		Fraction scale(1, int(rng()%5)+1);
		Q shift(Fraction(coord(),3),Fraction(coord(),3));
		auto map = [&](Q p) { return p*scale+shift; };
		auto from = [&](Z p) { return map(Q(p.x,p.y)); };
		Q a=from(za), b=from(zb), c=from(zc), d=from(zd);
		auto want = segmentOracle(za,zb,zc,zd);
		for (Q& p : want) p=map(p);
		auto got = segInter(a,b,c,d);
		sort(all(got)); sort(all(want)); assert(got==want);
		for (Q p : got) assert(onSegment(a,b,p) && onSegment(c,d,p));
		if (!(a==b)) {
			Q p = lineProj(a,b,c), r = lineProj(a,b,c,true);
			assert(a.cross(b,p)==0 && (p-c).dot(b-a)==0);
			assert((c+r)/2==p && (r-a).dist2()==(c-a).dist2());
			assert(polygon::linearTransformation(a,b,c,d,a)==c);
			assert(polygon::linearTransformation(a,b,c,d,b)==d);
		}
		if (a.cross(b,c)!=0) {
			Q center=polygon::ccCenter(a,b,c);
			assert((center-a).dist2()==(center-b).dist2());
			assert((center-a).dist2()==(center-c).dist2());
			assert(polygon::polygonCenter(vector<Q>{a,b,c})==(a+b+c)/3);
		}
		if (!(a==b) && !(c==d)) {
			auto p=lineInter(a,b,c,d);
			ll den=(zb-za).cross(zd-zc);
			if (den) {
				assert(p.first==1);
				assert(p.second==a+(b-a)*Fraction((zc-za).cross(zd-zc),den));
				assert(a.cross(b,p.second)==0 && c.cross(d,p.second)==0);
			} else assert(p.first==(za.cross(zb,zc)==0 ? -1 : 0));
		}
	}
	rep(it,0,2000) {
		vector<Q> points;
		rep(i,0,12) points.emplace_back(Fraction(coord(),3),Fraction(coord(),3));
		auto h=hull::convexHull(points);
		if (sz(h)<3) continue;
		for (Q p : points) {
			assert(hull::inHull(h,p,false));
			assert(polygon::inPolygon(h,p,false));
		}
		rep(i,0,sz(h)) assert(h[i].cross(h[(i+1)%sz(h)],h[(i+2)%sz(h)])>0);
		Fraction diameter=0;
		for (Q p : points) for (Q q : points)
			diameter=max(diameter,(p-q).dist2());
		auto ends=hull::hullDiameter(h);
		assert((ends[0]-ends[1]).dist2()==diameter);
		Q p(Fraction(coord(),3),Fraction(coord(),3));
		bool inside=true, strict=true;
		rep(i,0,sz(h)) {
			Fraction c=h[i].cross(h[(i+1)%sz(h)],p);
			inside &= c>=0; strict &= c>0;
		}
		assert(hull::inHull(h,p,false)==inside);
		assert(hull::inHull(h,p,true)==strict);
		assert(polygon::inPolygon(h,p,false)==inside);
		assert(polygon::inPolygon(h,p,true)==strict);
		Q dir(Fraction(1,2),Fraction(2,3));
		int e=extrVertex(h,dir);
		for (Q q : h) assert(h[e].dot(dir)>=q.dot(dir));
		auto left=polygon::polygonCut(h,p,p+dir);
		auto right=polygon::polygonCut(h,p+dir,p);
		assert(area(left)+area(right)==area(h));
		for (Q q : left) assert(p.cross(p+dir,q)<=0);
		for (Q q : right) assert(p.cross(p+dir,q)>=0);
		vector<halfplane::RayT<Fraction>> rays;
		rep(i,0,sz(h)) rays.push_back({h[i],h[(i+1)%sz(h)]-h[i]});
		shuffle(all(rays),rng);
		auto hp=halfplane::halfPlaneIsect(rays);
		assert(area(hp)==area(h));
		if (area(right)!=0) {
			rays.push_back({p,dir});
			hp=halfplane::halfPlaneIsect(rays);
			assert(area(hp)==area(right));
		}
		vector<Q> hits;
		rep(i,0,sz(h)) {
			auto cross=lineInter(p,p+dir,h[i],h[(i+1)%sz(h)]);
			if (cross.first==1 && onSegment(h[i],h[(i+1)%sz(h)],cross.second))
				hits.push_back(cross.second);
		}
		auto hit=lineHull(p,p+dir,h);
		assert((hit[0]!=-1)==!hits.empty());
	}
	rep(it,0,2000) {
		Fraction x(coord(),3), y(coord(),3);
		Fraction w(int(rng()%7)+1,3), h(int(rng()%7)+1,3);
		auto a=rect(0,0,2,3), b=rect(x,y,w,h);
		Fraction dx=max(Fraction(0),min(Fraction(2),x+w)-max(Fraction(0),x));
		Fraction dy=max(Fraction(0),min(Fraction(3),y+h)-max(Fraction(0),y));
		vector<vector<Q>> polys={a,b};
		assert(polygon::polyUnion(polys)==6+w*h-dx*dy);
		polys.push_back(a); assert(polygon::polyUnion(polys)==6+w*h-dx*dy);
	}
	using Q3=Point3D<Fraction>;
	vector<Q3> tetra={Q3(),Q3(Fraction(1,2),0,0),
		Q3(0,Fraction(2,3),0),Q3(0,0,Fraction(3,4))};
	struct Face { int a,b,c; };
	vector<Face> faces={{0,2,1},{0,1,3},{0,3,2},{1,2,3}};
	assert(signedPolyVolume(tetra,faces)==Fraction(1,24));
	Q3 cross=tetra[1].cross(tetra[2]);
	assert(cross.dot(tetra[1])==0 && cross.dot(tetra[2])==0);
	assert(abs(Q(3,4).dist()-5)<1e-12);
	assert(abs(Q(0,1).angle()-acos(-1)/2)<1e-12);
	assert(abs(Q3(0,1,0).phi()-acos(-1)/2)<1e-12);
	assert(abs(Q3(0,0,1).theta())<1e-12);
	cout << "Fraction geometry tests passed (24000 random cases)\n";
}
