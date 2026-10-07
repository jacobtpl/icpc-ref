#include "../utilities/template.h"

// from content/contest/template.cpp
template<class A, class B> bool ckmin(A& a, const B& b) { return b < a ? a = b, 1 : 0; }

#include "../../content/graph/BipolarOrientation.h"

mt19937 rng(55555);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

bool connectedWithout(const vector<vi>& a, int del) {
	int n = sz(a), start = del == 0 ? 1 : 0;
	if (n - (del >= 0) <= 0) return 1;
	vi seen(n), q = {start};
	seen[start] = 1;
	rep(i,0,sz(q)) for (int y : a[q[i]]) if (y != del && !seen[y])
		seen[y] = 1, q.push_back(y);
	return sz(q) == n - (del >= 0);
}
bool biconnected(const vector<vi>& a) {
	rep(del,-1,sz(a)) if (!connectedWithout(a, del)) return 0;
	return 1;
}

// st-numbering check: s first, t last, every other vertex has a neighbour
// on both sides (so orienting edges low -> high gives a bipolar orientation)
void check(const vector<vi>& a, int s, int t) {
	int n = sz(a);
	vector<vi> copy = a;
	vi ord = bipolarOrient(copy, s, t);
	assert(copy == a);
	assert(sz(ord) == n && ord[0] == s && ord[n-1] == t);
	vi pos(n, -1);
	rep(i,0,n) {
		assert(0 <= ord[i] && ord[i] < n && pos[ord[i]] == -1);
		pos[ord[i]] = i;
	}
	rep(v,0,n) if (v != s && v != t) {
		bool lo = 0, hi = 0;
		for (int x : a[v]) lo |= pos[x] < pos[v], hi |= pos[x] > pos[v];
		if (!lo || !hi) {
			cerr << "bad order n=" << n << " s=" << s << " t=" << t << " v=" << v << endl;
			rep(u,0,n) for (int x : a[u]) if (u <= x) cerr << u << ' ' << x << endl;
			for (int x : ord) cerr << x << ' ';
			cerr << endl;
			abort();
		}
	}
}

int main() {
	// K2, with and without a doubled edge
	check({{1}, {0}}, 0, 1);
	check({{1}, {0}}, 1, 0);
	check({{1, 1}, {0, 0}}, 0, 1);
	// small random biconnected graphs, every (s, t)
	int tested = 0;
	rep(it,0,200000) {
		int n = ri(3, 8), m = ri(n, min(n * (n-1) / 2 + 2, 2 * n));
		vector<vi> a(n);
		rep(i,0,m) {
			int u = ri(0, n-1), v = ri(0, n-2);
			if (v >= u) v++;
			if (ri(0, 9) && count(all(a[u]), v)) continue; // few multi-edges
			a[u].push_back(v);
			a[v].push_back(u);
		}
		if (!biconnected(a)) continue;
		tested++;
		rep(s,0,n) rep(t,0,n) if (s != t) check(a, s, t);
	}
	assert(tested > 20000);
	// larger: Hamiltonian cycle / ear-like structures plus random chords
	rep(it,0,2000) {
		int n = ri(3, 300), extra = ri(0, 2) ? ri(0, 3) : ri(0, 2 * n);
		vi perm(n);
		iota(all(perm), 0);
		shuffle(all(perm), rng);
		vector<vi> a(n);
		auto add = [&](int u, int v) { a[u].push_back(v), a[v].push_back(u); };
		rep(i,0,n) add(perm[i], perm[(i+1) % n]);
		rep(i,0,extra) {
			int u = ri(0, n-1), v = ri(0, n-2);
			if (v >= u) v++;
			add(u, v);
		}
		for (auto& v : a) shuffle(all(v), rng);
		rep(k,0,5) {
			int s = ri(0, n-1), t = ri(0, n-2);
			if (t >= s) t++;
			check(a, s, t);
		}
	}
	cout << "Tests passed!" << endl;
}
