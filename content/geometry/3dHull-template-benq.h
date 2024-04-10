/**
 * Author: Benq
 * Description: 3d hull helpers
 */

using T = double;
using P3 = Point3D<double>;
using vb = vector<bool>;

mt19937 rng;
P3 cross(const P3& a, const P3& b, const P3& c) {
	return (b-a).cross(c-a); }
P3 perp(const P3& a, const P3& b, const P3& c) {
	return cross(a,b,c).unit(); }
bool isMult(const P3& a, const P3& b) { // for long longs
	P3 c = a.cross(b);
	return (c.x== 0 && c.y==0 && c.z==0);
}
bool collinear(const P3& a, const P3& b, const P3& c) { 
	return isMult(b-a,c-a); }
T DC(const P3&a,const P3&b,const P3&c,const P3&p) { 
	return cross(a,b,c).dot(p-a); }
bool coplanar(const P3&a,const P3&b,const P3&c,const P3&p) { 
	return DC(a,b,c,p) == 0; }
bool above(const P3&a,const P3&b,const P3&c,const P3&p) { 
	return DC(a,b,c,p) > 0; } // is p strictly above plane
void prep(vector<P3>& p) { // rearrange points such that
	shuffle(all(p),rng); // first four are not coplanar
	int dim = 1; 
	rep(i,1,sz(p)) 
		if (dim == 1) {
			if (p[0] != p[i]) swap(p[1],p[i]), ++dim;
		} else if (dim == 2) {
			if (!collinear(p[0],p[1],p[i])) 
				swap(p[2],p[i]), ++dim;
		} else if (dim == 3) {
			if (!coplanar(p[0],p[1],p[2],p[i]))
				swap(p[3],p[i]), ++dim;
		}
	assert(dim == 4);
}
