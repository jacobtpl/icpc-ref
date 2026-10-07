#include "../utilities/template.h"

#include "../../content/graph/MaximumClique.h"
namespace maximal {
#include "../../content/graph/MaximalCliques.h"
}

struct timeit {
	decltype(chrono::high_resolution_clock::now()) begin;
	const string label;
	timeit(string label = "???") : label(label) { begin = chrono::high_resolution_clock::now(); }
	~timeit() {
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
		cerr << duration << "ms elapsed [" << label << "]" << endl;
	}
};


bool isClique(vb& ed, vi c) {
	sort(all(c));
	rep(i,0,sz(c)) rep(j,0,i)
		if (c[i] == c[j] || !ed[c[i]][c[j]]) return false;
	for (int x : c) if (x < 0 || x >= sz(ed)) return false;
	return true;
}

// Exhaustive bitmask oracle for n <= 16; checks the returned set itself.
void testBrute() {
	rep(it,0,36000) {
		int n = rand() % (it < 30000 ? 7 : 16) + 1, p = rand() % 101;
		vb ed(n); vi adj(n);
		rep(i,0,n) rep(j,0,i) if (rand() % 100 < p)
			ed[i][j] = ed[j][i] = 1, adj[i] |= 1 << j, adj[j] |= 1 << i;
		int best = 0;
		rep(m,1,1 << n) {
			bool ok = 1;
			rep(i,0,n) if (m >> i & 1) ok &= (m & ~adj[i] & ~(1 << i)) == 0;
			if (ok) best = max(best, __builtin_popcount(m));
		}
		vi c = Maxclique(ed).maxClique();
		assert(isClique(ed, c) && sz(c) == best);
	}
}

// Structured graphs up to the full bitset<200> width with known answers.
void testBig() {
	rep(it,0,300) {
		int n = it < 20 ? 200 : rand() % 200 + 1;
		// complete k-partite graph: max clique = number of parts
		int k = rand() % n + 1;
		vi col(n);
		rep(i,0,n) col[i] = i < k ? i : rand() % k;
		random_shuffle(all(col));
		vb ed(n);
		rep(i,0,n) rep(j,0,i) if (col[i] != col[j]) ed[i][j] = ed[j][i] = 1;
		vi c = Maxclique(ed).maxClique();
		assert(isClique(ed, c) && sz(c) == k);
		// disjoint cliques of random sizes: max clique = largest block
		vi cnt(k);
		rep(i,0,n) cnt[col[i]]++;
		rep(i,0,n) rep(j,0,i) ed[i][j] = ed[j][i] = col[i] == col[j];
		c = Maxclique(ed).maxClique();
		assert(isClique(ed, c) && sz(c) == *max_element(all(cnt)));
	}
}

int main() {
    srand(6);
    testBrute();
    testBig();
#ifdef TEST_EMPTY
    // Opt-in: n = 0 reads r[0] in init() and S[1] in expand() (out of
    // bounds; caught by -D_GLIBCXX_DEBUG / ASan).
    assert(Maxclique(vb()).maxClique().empty());
#endif
    rep(it, 0, 100000) {
        int n =(rand()%32)+1;
        vb ed(n);
        vector<maximal::B> ed2(n);
        int p =rand()%100;
        rep(i, 0, n) rep(j, 0, i) {
            ed[i][j] = (rand() % 100) < p;
            ed[j][i] = ed[i][j];
            ed2[i][j] = ed[i][j];
            ed2[j][i] = ed[j][i];
        }
        Maxclique clique2(ed);
        int mx = 0;
        maximal::cliques(ed2, [&](auto x){mx = max(mx, int(x.count()));});
        assert(mx == sz(clique2.maxClique()));
    }
    cout<<"Tests passed!"<<endl;
}

