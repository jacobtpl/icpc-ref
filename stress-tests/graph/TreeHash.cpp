#include "../utilities/template.h"
#include "../utilities/genTree.h"

#include "../../content/graph/TreeHash.h"

// Oracle: canonical parenthesis string of the subtree of v rooted at p.
string canon(vector<vi>& g, int v, int p) {
	vector<string> c;
	for (int u : g[v]) if (u != p) c.push_back(canon(g, u, v));
	sort(all(c));
	string s = "(";
	for (auto& x : c) s += x;
	return s + ")";
}
string unrootedCanon(vector<vi>& g) {
	string best;
	rep(v,0,sz(g)) { string s = canon(g, v, -1); if (!v || s < best) best = s; }
	return best;
}
// Oracle 2: brute force over all vertex bijections (n <= 7).
bool bruteIso(vector<vi>& a, vector<vi>& b) {
	int n = sz(a);
	if (n != sz(b)) return false;
	set<pii> eb;
	rep(v,0,n) for (int u : b[v]) eb.insert({v, u});
	vi perm(n); iota(all(perm), 0);
	do {
		bool ok = true;
		rep(v,0,n) for (int u : a[v]) if (!eb.count({perm[v], perm[u]})) { ok = false; break; }
		if (ok) return true;
	} while (next_permutation(all(perm)));
	return false;
}

vector<vi> fromEdges(int n, const vector<pii>& e) {
	vector<vi> g(n);
	for (auto [a, b] : e) g[a].push_back(b), g[b].push_back(a);
	return g;
}
vector<vi> relabel(vector<vi>& g) {
	int n = sz(g);
	mt19937 rng(rand());
	vi p(n); iota(all(p), 0); shuffle(all(p), rng);
	vector<vi> h(n);
	rep(v,0,n) for (int u : g[v]) h[p[v]].push_back(p[u]);
	for (auto& a : h) shuffle(all(a), rng);
	return h;
}
vector<vi> randTree(int n) {
	if (n == 1) return vector<vi>(1);
	int t = rand() % 6;
	vector<pii> e;
	if (t == 0) rep(i,1,n) e.push_back({i - 1, i}); // path
	else if (t == 1) rep(i,1,n) e.push_back({0, i}); // star
	else if (t == 2) rep(i,1,n) e.push_back({(i - 1) / 2, i}); // binary
	else if (t == 3) rep(i,1,n) e.push_back({rand() % i, i}); // random recursive
	else if (t == 4) rep(i,1,n) e.push_back({max(0, i - 1 - rand() % 3), i}); // low branching
	else e = genRandomTree(n); // uniform (Prufer)
	auto g = fromEdges(n, e);
	return relabel(g);
}

void checkSiz(vector<vi>& g, int r) {
	treeId(g, r);
	vi s = siz;
	// recompute sizes with the oracle (string length / 2 = subtree size)
	function<void(int,int)> dfs = [&](int v, int p) {
		assert(s[v] == sz(canon(g, v, p)) / 2);
		for (int u : g[v]) if (u != p) dfs(u, v);
	};
	dfs(r, -1);
}

void testRootedPool(int maxN, int trees, int iters) {
	rep(it,0,iters) {
		tid.clear();
		vector<vector<vi>> gs;
		vector<vi> ids;
		vector<vector<string>> cs; // canonical string of subtree of v (rooted at root)
		vi roots;
		rep(t,0,trees) {
			int n = rand() % maxN + 1;
			gs.push_back(randTree(n));
			roots.push_back(rand() % n);
		}
		if (trees > 1 && rand() % 2) { gs[1] = relabel(gs[0]); roots[1] = roots[0]; }
		rep(t,0,trees) {
			ids.push_back(treeId(gs[t], roots[t]));
			vector<string> c(sz(gs[t]));
			function<void(int,int)> dfs = [&](int v, int p) {
				c[v] = canon(gs[t], v, p);
				for (int u : gs[t][v]) if (u != p) dfs(u, v);
			};
			dfs(roots[t], -1);
			cs.push_back(c);
			// every id must be < sz(tid) and ids must be dense 0..sz(tid)-1 overall
			for (int x : ids[t]) assert(0 <= x && x < sz(tid));
		}
		set<int> used;
		rep(a,0,trees) rep(v,0,sz(gs[a])) {
			used.insert(ids[a][v]);
			rep(b,0,trees) rep(u,0,sz(gs[b]))
				assert((ids[a][v] == ids[b][u]) == (cs[a][v] == cs[b][u]));
		}
		assert(sz(used) == sz(tid));
	}
}

void testUnrooted(int maxN, int iters, bool brute) {
	rep(it,0,iters) {
		if (rand() % 3 == 0) tid.clear();
		int n = rand() % maxN + 1;
		auto g1 = randTree(n);
		auto g2 = rand() % 2 ? relabel(g1) : randTree(n);
		int u1 = unrootedId(g1), u2 = unrootedId(g2);
		bool iso = unrootedCanon(g1) == unrootedCanon(g2);
		assert((u1 == u2) == iso);
		if (brute) assert(bruteIso(g1, g2) == iso);
		assert(unrootedId(g1) == u1); // idempotent
		auto g3 = relabel(g1);
		assert(unrootedId(g3) == u1);
		// unrooted id is a rooted id of some root, and its tid entry is shared
		bool found = false;
		rep(v,0,n) if (treeId(g1, v)[v] == u1) found = true;
		assert(found);
		assert(sz(tid) > u1);
	}
}

// Exhaustive: all labelled trees on n vertices via Prufer codes must be
// partitioned by unrootedId exactly as by the string oracle.
void testExhaustive(int n) {
	tid.clear();
	map<string, int> cls;
	int len = n - 2, total = 1;
	rep(i,0,len) total *= n;
	rep(code,0,total) {
		vi pc(len); int c = code;
		rep(i,0,len) pc[i] = c % n + 1, c /= n;
		auto e = pruferCodeToTree(pc);
		for (auto& p : e) p.first--, p.second--;
		auto g = fromEdges(n, e);
		int id = unrootedId(g);
		string s = unrootedCanon(g);
		auto it = cls.find(s);
		if (it == cls.end()) cls[s] = id;
		else assert(it->second == id);
	}
	set<int> distinct;
	for (auto& [s, id] : cls) distinct.insert(id);
	assert(sz(distinct) == sz(cls));
	static const int unlabelled[] = {0, 1, 1, 1, 2, 3, 6, 11, 23};
	assert(sz(cls) == unlabelled[n]);
}

void testEdgeCases() {
	tid.clear();
	vector<vi> one(1);
	assert(sz(treeId(one, 0)) == 1 && siz == vi{1});
	int i1 = unrootedId(one);
	assert(i1 == treeId(one, 0)[0]);
	vector<vi> two = fromEdges(2, {{0, 1}});
	vi id2 = treeId(two, 1);
	assert(id2[1] != id2[0] && id2[0] == i1 && siz[1] == 2 && siz[0] == 1);
	assert(unrootedId(two) == treeId(two, 0)[0] && unrootedId(two) != i1);
	// a leaf anywhere has the same id as the single vertex tree
	auto g = randTree(50);
	vi id = treeId(g, 7);
	rep(v,0,50) if (v != 7 && sz(g[v]) == 1) assert(id[v] == i1);
	// tid.clear() resets ids to start from 0 again
	tid.clear();
	assert(treeId(one, 0)[0] == 0 && sz(tid) == 1);
	// ids survive across calls: same tree twice gives identical ids
	auto h = randTree(200);
	assert(treeId(h, 3) == treeId(h, 3));
	int before = sz(tid);
	treeId(h, 3);
	assert(sz(tid) == before);
}

void testLarge() {
	tid.clear();
	int n = 1000000; // no recursion: deep path must not overflow the stack
	vector<pii> e; rep(i,1,n) e.push_back({i - 1, i});
	auto g = fromEdges(n, e);
	vi id = treeId(g, 0);
	rep(v,0,n) assert(id[v] == n - 1 - v && siz[v] == n - v);
	assert(sz(tid) == n);
	tid.clear();
	n = 200000;
	e.clear(); rep(i,1,n) e.push_back({i - 1, i});
	auto p = fromEdges(n, e);
	auto p2 = relabel(p);
	assert(unrootedId(p) == unrootedId(p2));
	e.clear(); rep(i,1,n) e.push_back({0, i});
	auto s = fromEdges(n, e);
	auto s2 = relabel(s);
	assert(unrootedId(s) == unrootedId(s2) && unrootedId(s) != unrootedId(p));
	e.clear(); rep(i,1,n) e.push_back({rand() % i, i});
	auto r = fromEdges(n, e);
	auto r2 = relabel(r);
	assert(unrootedId(r) == unrootedId(r2) && unrootedId(r) != unrootedId(p));
}

int main() {
	srand(2);
	testEdgeCases();
	testUnrooted(7, 300, true);
	testUnrooted(40, 300, false);
	testRootedPool(12, 4, 300);
	testRootedPool(40, 3, 100);
	rep(i,0,100) { auto g = randTree(rand() % 30 + 1); checkSiz(g, rand() % sz(g)); }
	rep(n,3,8) testExhaustive(n);
	testLarge();
	cout << "Tests passed!" << endl;
}
