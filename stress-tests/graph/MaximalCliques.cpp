#include "../utilities/template.h"

#include "../../content/graph/MaximalCliques.h"

template<class F>
void fastCliques(vector<B>& eds, F f) {
	B R{}, P = ~B(), X{};
	vi deg(sz(eds));
	rep(i,0,sz(eds)) deg[i] = sz(eds[i]);
	rep(j,0,sz(eds)) {
		int i = (int)(min_element(all(deg)) - deg.begin());
		R[i] = 1;
		rec(eds, R, P & eds[i], X & eds[i], f);
		R[i] = P[i] = 0; X[i] = 1;
		rep(k,0,sz(eds)) if (eds[i][k]) deg[k]--;
		deg[i] = 1000000;
	}
}

int main1() {
	rep(n,1,11) rep(m,0,200) {
		vector<B> ed(n);
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n;
			if (a == b) continue;
			ed[a][b] = 1;
			ed[b][a] = 1;
		}
		unordered_set<B> cl;
		int co = 0;
		cliques(ed, [&](B x) {
			co++;
			cl.insert(x);
		});
		assert(sz(cl) == co); // no duplicates
		auto isClique = [&](B c) {
			rep(i,0,n) if (c[i])
			rep(j,i+1,n) if (c[j]) {
				if (!ed[i][j]) return false; // not a clique
			}
			rep(i,0,n) if (!c[i]) {
				bool all = true;
				rep(j,0,n) if (c[j]) all &= ed[i][j];
				if (all) return false; // not maximal
			}
			return true;
		};
		for(auto &c: cl) {
			assert(isClique(c));
		}

		int realCo = 0;
		rep(bi,0,(1 << n)) {
			B c{};
			rep(i,0,n) c[i] = !!(bi & (1 << i));
			if (isClique(c)) realCo++;
		}
		assert(co == realCo);
	}
	// n up to 16, all densities, set-equality against bitmask brute force
	srand(7);
	rep(it,0,3000) {
		int n = rand() % (it < 2700 ? 12 : 16) + 1, p = rand() % 101;
		vector<B> ed(n); vi adj(n);
		rep(i,0,n) rep(j,0,i) if (rand() % 100 < p)
			ed[i][j] = ed[j][i] = 1, adj[i] |= 1 << j, adj[j] |= 1 << i;
		set<unsigned long long> got, want;
		int co = 0;
		cliques(ed, [&](B x) { co++; got.insert(x.to_ullong()); });
		assert(co == sz(got));
		rep(m,1,1 << n) {
			int common = (1 << n) - 1;
			rep(i,0,n) if (m >> i & 1) common &= adj[i];
			// clique iff every member is adjacent to all others;
			// maximal iff no outside vertex is adjacent to all members
			bool ok = 1;
			rep(i,0,n) if (m >> i & 1) ok &= (m & ~adj[i] & ~(1 << i)) == 0;
			if (ok && !common) want.insert(m);
		}
		assert(got == want);
	}
	// Full bitset<128> width: disjoint cliques (one maximal clique per
	// block) and complete multipartite graphs (product of part sizes).
	rep(it,0,200) {
		int n = it < 10 ? 128 : rand() % 128 + 1, k = rand() % n + 1;
		vi col(n), cnt(k);
		rep(i,0,n) col[i] = i < k ? i : rand() % k;
		random_shuffle(all(col));
		rep(i,0,n) cnt[col[i]]++;
		vector<B> ed(n);
		rep(i,0,n) rep(j,0,i) ed[i][j] = ed[j][i] = col[i] == col[j];
		int co = 0;
		cliques(ed, [&](B x) {
			co++;
			int c = col[(int)x._Find_first()];
			assert((int)x.count() == cnt[c]);
			rep(i,0,n) assert(x[i] == (col[i] == c));
		});
		assert(co == k);
		double prod = 1;
		for (int c : cnt) prod *= c;
		if (prod > 2e5) continue;
		rep(i,0,n) rep(j,0,i) ed[i][j] = ed[j][i] = col[i] != col[j];
		co = 0;
		cliques(ed, [&](B x) { co++; assert((int)x.count() == k); });
		assert(co == (int)prod);
	}
#ifdef TEST_EMPTY
	{ // Opt-in: n = 0 reads eds[0] (segfault at -O2)
		vector<B> ed;
		int co = 0;
		cliques(ed, [&](B x) { co++; assert(x.none()); });
		assert(co <= 1);
	}
#endif
	cout<<"Tests passed!"<<endl;
	return 0;
}

int main2() {
	rep(it,0,20) {
		const int n = 128, m = 4000;
		vector<B> ed(n);
		rep(i,0,m) {
			int a = rand() % n, b = rand() % n;
			if (a == b) continue;
			ed[a][b] = 1;
			ed[b][a] = 1;
		}
		int co = 0, sum = 0;
		cliques(ed, [&](B x) { co++; sum += (int)x.count(); });
		cout << co << ' ' << (double)sum / co << endl;
	}
	return 0;
}

#ifndef target
#define target main1
#endif
int main() { target(); }
