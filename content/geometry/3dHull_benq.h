/**
 * Author: Benq
 * Description: Computes all faces of the 3-dimension hull of a point set.
 *  *No four points must be coplanar*, or else random results will be returned.
 *  All faces will point outwards.
 * Time: O(n^2, n \log n)
 * Status: tested on SPOJ CH3D
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

using F = array<int,3>; // face
vector<F> hull3d(vector<P3>& p) {
	// s.t. first four points form tetra
	prep(p); int N = sz(p); vector<F> hull; // triangle for each face
	auto ad = [&](int a, int b, int c) { hull.pb({a,b,c}); }; 
	// +new face to hull
	ad(0,1,2), ad(0,2,1); // initialize hull as first 3 points
	vector<vb> in(N,vb(N)); // is zero before each iteration
	rep(i,3,N) { // incremental construction
		vector<F> def, HULL; swap(hull,HULL); 
		// HULL now contains old hull
		auto ins = [&](int a, int b, int c) {
			if (in[b][a]) in[b][a] = 0; // kill reverse face
			else in[a][b] = 1, ad(a,b,c);
		};
		for (auto &f:HULL) {
			if (above(p[f[0]],p[f[1]],p[f[2]],p[i])) 
				rep(j,0,3) ins(f[j],f[(j+1)%3],i); 
				// recalc all faces s.t. point is above face
			else def.pb(f); 
		}
		for (auto &t:hull) if (in[t[0]][t[1]]) // edge exposed, 
			in[t[0]][t[1]] = 0, def.pb(t); // add a new face
		swap(hull,def);
	}
	return hull;
}
vector<F> hull3dFast(vector<P3>& p) {
	prep(p); int N = sz(p); vector<F> hull; 
	vb active; // whether face is active
	vector<vi> rvis; // points visible from each face
	vector<array<pi,3>> other; // other face adjacent to each edge of face
	vector<vi> vis(N); // faces visible from each point
	auto ad = [&](int a, int b, int c) { 
		hull.pb({a,b,c}); active.pb(1); rvis.emplace_back(); other.emplace_back(); };
	auto ae = [&](int a, int b) { vis[b].pb(a), rvis[a].pb(b); };
	auto abv = [&](int a, int b) {
		F f=hull[a]; return above(p[f[0]],p[f[1]],p[f[2]],p[b]);};
	auto edge = [&](pi e) -> pi { 
		return {hull[e.first][e.second],hull[e.first][(e.second+1)%3]}; };
	auto glue = [&](pi a, pi b) { // link two faces by an edge
		pi x = edge(a); assert(edge(b) == mp(x.second,x.first));
		other[a.first][a.second] = b, other[b.first][b.second] = a;
	}; // ensure face 0 is removed when i=3
	ad(0,1,2), ad(0,2,1); if (abv(1,3)) swap(p[1],p[2]); 
	rep(i,0,3) glue({0,i},{1,2-i});
	rep(i,3,N) ae(abv(1,i),i); // coplanar points go in rvis[0]
	vi label(N,-1);
	rep(i,3,N) { // incremental construction
		vi rem; for(auto &t:vis[i]) if (active[t]) active[t]=0, rem.pb(t);
		if (!sz(rem)) continue; // hull unchanged
		int st = -1; 
		for(auto &r:rem) rep(j,0,3) {
			int o = other[r][j].first;
			if (active[o]) { // create new face!
				int a,b; tie(a,b) = edge({r,j}); ad(a,b,i); st = a;
				int cur = sz(rvis)-1; label[a] = cur; 
				vi tmp; set_union(all(rvis[r]),all(rvis[o]),
									back_inserter(tmp)); 
				// merge sorted vectors ignoring duplicates
				for(auto &x:tmp) if (abv(cur,x)) ae(cur,x);
				/// if no rounding errors then guaranteed that only x>i matters
				glue({cur,0},other[r][j]); // glue old w/ new face
			}
		}
		for (int x = st, y; ; x = y) { // glue new faces together
			int X = label[x]; glue({X,1},{label[y=hull[X][1]],2});
			if (y == st) break;
		}
	}
	vector<F> ans; rep(i,0,sz(hull)) if (active[i]) ans.pb(hull[i]);
	return ans;
}