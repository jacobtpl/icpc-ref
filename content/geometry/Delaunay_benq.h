/**
 * Author: Benq
 * Date: Unknown
 * Description: Computes the Delaunay triangulation of a set of points.
 *  Each circumcircle contains none of the input points.
 * Time: O(hull3d)
 * Status: https://dmoj.ca/problem/cco08p6
 */

using P = Point<ll>;
vector<array<P,3>> triHull(vector<P> p) {
	vector<P3> p3; vector<array<P,3>> res; for (auto &x:p) p3.pb(P3{x.x,x.y,x.dist2()});
	bool ok = 0; for (auto &t:p3) ok |= !coplanar(p3[0],p3[1],p3[2],t);
	if (!ok) { // all points concyclic
		sort(1+all(p),[&p](P a, P b) { 
			return (a-p.front()).cross(b-p.front())>0; });
		rep(i,1,sz(p)-1) res.pb({p.front(),p[i],p[i+1]});
	} else {
		#define nor(z) P(p3[z].x,p3[z].y)
		for(auto &t:hull3dFast(p3)) 
			if (cross(p3[t[0]],p3[t[1]],p3[t[2]]).dot(P3{0,0,1}) < 0)
				res.pb({nor(t[0]),nor(t[2]),nor(t[1])});
	}
	return res;
}