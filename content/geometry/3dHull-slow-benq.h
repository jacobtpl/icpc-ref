/**
 * Author: Benq
 * Description: Computes all faces of the 3-dimension hull of a point set.
 *  *No four points must be coplanar*, or else random results will be returned.
 *  All faces will point outwards.
 * Time: O(n^2)
 * Status: tested on SPOJ CH3D
 */

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
