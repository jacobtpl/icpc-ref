#include "../utilities/template.h"

// content/graph/DynamicTree.h is not a pasteable header: it is a complete
// solution program (a top tree over a link-cut tree, with path and subtree
// assign/add and min/max/sum queries under rerooting and reparenting) that
// depends on a competitive programming template which is not in this repo
// (MX, vpi, setIO, re, ps, pr) and that has its own main() reading stdin.
// This test supplies the missing pieces, renames main() and feeds it through
// string streams, comparing against a brute force.
//
// Input format of the program:
//   N M, N-1 edges, N values, root, then M operations
//   0 x y   assign y to the subtree of x      5 x y   add y to the subtree of x
//   1 x     make x the root                   9 x y   make y the parent of x
//   2 x y z assign z to the path x..y         6 x y z add z to the path x..y
//   3 x / 4 x / 11 x   min / max / sum of the subtree of x
//   7 x y / 8 x y / 10 x y   min / max / sum of the path x..y

mt19937 rng(987654321);
int rnd(int lo, int hi) { return (int)(rng() % (unsigned)(hi - lo + 1)) + lo; }

struct Brute {
	int n, root;
	vector<set<int>> adj;
	vi val, par;
	void dfs(int v, int p) {
		par[v] = p;
		for (int w : adj[v]) if (w != p) dfs(w, v);
	}
	void subtree(int v, vi& out) {
		out.push_back(v);
		for (int w : adj[v]) if (w != par[v]) subtree(w, out);
	}
	vi nodes(int k, int x, int y) {
		vi r;
		if (k == 0 || k == 5 || k == 3 || k == 4 || k == 11) {
			dfs(root, 0);
			subtree(x, r);
		} else {
			dfs(x, 0);
			for (; y; y = par[y]) r.push_back(y);
		}
		return r;
	}
	// returns true and sets res if the operation prints something
	bool apply(int k, int x, int y, int z, ll& res) {
		if (k == 1) { root = x; return false; }
		if (k == 9) {
			dfs(root, 0);
			for (int v = y; v; v = par[v]) if (v == x) return false;
			adj[x].erase(par[x]), adj[par[x]].erase(x);
			adj[x].insert(y), adj[y].insert(x);
			return false;
		}
		vi vs = nodes(k, x, y);
		if (k == 0) for (int v : vs) val[v] = y;
		else if (k == 5) for (int v : vs) val[v] += y;
		else if (k == 2) for (int v : vs) val[v] = z;
		else if (k == 6) for (int v : vs) val[v] += z;
		else {
			ll mn = LLONG_MAX, mx = LLONG_MIN, sum = 0;
			for (int v : vs) mn = min(mn, (ll)val[v]), mx = max(mx, (ll)val[v]), sum += val[v];
			res = k == 3 || k == 7 ? mn : k == 4 || k == 8 ? mx : sum;
			return true;
		}
		return false;
	}
};

const int MX = 200005;
typedef vector<pii> vpi;
stringstream IN;
void setIO() {}
template<class T> void re(T& x) { IN >> x; }
template<class A, class B> void re(pair<A, B>& p) { IN >> p.first >> p.second; }
template<class T> void re(vector<T>& v) { for (auto& x : v) re(x); }
template<class T, class... U> void re(T& t, U&... u) { re(t); re(u...); }
template<class... T> void ps(const T&...) {}
template<class T> void pr(const T&) {}

#undef sz
#undef all
#define main() dynamicTreeUnused = 0; void dynamicTreeMain()
#include "../../content/graph/DynamicTree.h"
#undef main
#undef f
#undef s

struct Op { int k, x, y, z; };

vector<ll> runHeader(int n, const vector<pii>& ed, const vi& val, int root, const vector<Op>& ops) {
	IN.str(""); IN.clear();
	IN << n << ' ' << ops.size() << '\n';
	for (auto& e : ed) IN << e.first << ' ' << e.second << '\n';
	rep(i,1,n+1) IN << val[i] << ' ';
	IN << root << '\n';
	for (auto& o : ops) {
		IN << o.k << ' ' << o.x;
		if (o.k != 1 && o.k != 3 && o.k != 4 && o.k != 11) IN << ' ' << o.y;
		if (o.k == 2 || o.k == 6) IN << ' ' << o.z;
		IN << '\n';
	}
	stringstream out;
	auto old = cout.rdbuf(out.rdbuf());
	dynamicTreeMain();
	cout.rdbuf(old);
	rep(i,1,n+1) delete LCT[i];
	vector<ll> res;
	ll x;
	while (out >> x) res.push_back(x);
	return res;
}

vector<pii> genTree(int n, int shape) {
	vector<pii> ed;
	vi perm(n + 1);
	iota(all(perm), 0);
	shuffle(perm.begin() + 1, perm.end(), rng);
	rep(i,2,n+1) {
		int p = shape == 0 ? rnd(1, i - 1) : shape == 1 ? i - 1 : shape == 2 ? 1 : max(1, i - rnd(1, 3));
		int a = perm[i], b = perm[p];
		if (rnd(0, 1)) swap(a, b);
		ed.push_back({a, b});
	}
	shuffle(all(ed), rng);
	return ed;
}

Op genOp(int n, int maxv, int mix) {
	static const int kinds[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	static const int queries[] = {3, 4, 7, 8, 10, 11};
	int k = kinds[rnd(0, 11)];
	if (mix == 1 && rnd(0, 1)) k = queries[rnd(0, 5)];
	if (mix == 2 && rnd(0, 1)) k = rnd(0, 1) ? 9 : 1;
	return {k, rnd(1, n), rnd(k == 9 || (k >= 6 && k != 11) ? 1 : -maxv, k == 9 || (k >= 6 && k != 11) ? n : maxv), rnd(-maxv, maxv)};
}

void testCorrect() {
	rep(it,0,40000) {
		int n = it % 7 == 0 ? rnd(1, 3) : rnd(1, 12);
		int m = rnd(1, 40), maxv = it % 3 ? 10 : 100000;
		auto ed = genTree(n, rnd(0, 3));
		Brute b;
		b.n = n, b.root = rnd(1, n);
		b.adj.assign(n + 1, {}), b.val.assign(n + 1, 0), b.par.assign(n + 1, 0);
		for (auto& e : ed) b.adj[e.first].insert(e.second), b.adj[e.second].insert(e.first);
		rep(i,1,n+1) b.val[i] = rnd(-maxv, maxv);
		vi val0 = b.val;
		int root0 = b.root;
		vector<Op> ops;
		vector<ll> expected;
		rep(i,0,m) {
			Op o = genOp(n, maxv, it % 3);
			// keep 2/6 well-formed: y is a node for path operations
			if (o.k == 2 || o.k == 6) o.y = rnd(1, n);
			ops.push_back(o);
			ll r;
			if (b.apply(o.k, o.x, o.y, o.z, r)) expected.push_back(r);
		}
		vector<ll> got = runHeader(n, ed, val0, root0, ops);
		if (got != expected) {
			cerr << "mismatch: n=" << n << " root=" << root0 << endl;
			for (auto& e : ed) cerr << e.first << ' ' << e.second << endl;
			rep(i,1,n+1) cerr << val0[i] << ' ';
			cerr << endl;
			for (auto& o : ops) cerr << o.k << ' ' << o.x << ' ' << o.y << ' ' << o.z << endl;
			cerr << "expected:"; for (ll x : expected) cerr << ' ' << x;
			cerr << "\ngot:"; for (ll x : got) cerr << ' ' << x;
			cerr << endl;
			assert(0);
		}
	}
}

// larger inputs: only checks that it terminates and prints the right number of answers
void testLarge(int n, int m, int shape, bool print) {
	auto ed = genTree(n, shape);
	vi val(n + 1);
	rep(i,1,n+1) val[i] = rnd(-1000, 1000);
	vector<Op> ops;
	size_t queries = 0;
	rep(i,0,m) {
		Op o = genOp(n, 1000, 0);
		if (o.k == 2 || o.k == 6) o.y = rnd(1, n);
		if (o.k == 5 || o.k == 6) o.y = o.k == 5 ? rnd(-1, 1) : o.y, o.z = rnd(-1, 1);
		ops.push_back(o);
		queries += o.k == 3 || o.k == 4 || o.k == 7 || o.k == 8 || o.k == 10 || o.k == 11;
	}
	auto t0 = chrono::steady_clock::now();
	auto res = runHeader(n, ed, val, 1, ops);
	double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
	assert(res.size() == queries);
	if (print) cerr << "n=" << n << " m=" << m << " shape=" << shape << ": " << el << "s" << endl;
}

int main(int argc, char**) {
	testCorrect();
	rep(shape,0,4) testLarge(20000, 20000, shape, argc > 1);
	if (argc > 1) rep(shape,0,4) testLarge(100000, 100000, shape, true);
	cout << "Tests passed!" << endl;
}
