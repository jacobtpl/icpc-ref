#include "../utilities/template.h"

#include "../../content/graph/MaximumClique.h"

// MaximumIndependentSet.h is a note: "find a max clique of the complement".
// Check that recipe against a bitmask brute force.
int main() {
	srand(8);
	rep(it,0,36000) {
		int n = rand() % (it < 30000 ? 7 : 16) + 1, p = rand() % 101;
		vi adj(n);
		vb comp(n);
		rep(i,0,n) rep(j,0,i) {
			if (rand() % 100 < p) adj[i] |= 1 << j, adj[j] |= 1 << i;
			else comp[i][j] = comp[j][i] = 1;
		}
		int best = 0;
		rep(m,1,1 << n) {
			bool ok = 1;
			rep(i,0,n) if (m >> i & 1) ok &= !(m & adj[i]);
			if (ok) best = max(best, __builtin_popcount(m));
		}
		vi c = Maxclique(comp).maxClique();
		assert(sz(c) == best);
		int m = 0;
		for (int x : c) assert(0 <= x && x < n && !(m >> x & 1)), m |= 1 << x;
		for (int x : c) assert(!(m & adj[x]));
	}
	cout<<"Tests passed!"<<endl;
}
