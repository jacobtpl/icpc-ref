/**
 * Author: Benq
 * Description: Computes all faces of the 3-dimension hull of a point set.
 * Time: O(n \log n)
 * Status: tested on SPOJ CH3D
 */

using F = array<int,3>; // face
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
