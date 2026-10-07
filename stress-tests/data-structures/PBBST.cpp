#include "../utilities/template.h"

#include "../../content/data-structures/PBBST.h"

// Usage contract exercised here (the header documents none): nodes created
// since the last step() are owned by exactly one handle and are consumed by
// split / operator+. step() freezes every existing node, so all handles that
// exist at that moment become immutable versions that may be reused freely;
// the handle step() returns is a fresh mutable copy (its root belongs to the
// new epoch), so it is again consumed by the next split / operator+.

mt19937 rng(987654321);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

void collect(int n, vi& out) {
	if (n == -1) return;
	collect(PAVL::N[n].c[0], out);
	out.push_back(PAVL::N[n].val);
	collect(PAVL::N[n].c[1], out);
}
// checks AVL invariant + cached size/height, returns height
int check(int n) {
	if (n == -1) return 0;
	auto& x = PAVL::N[n];
	int a = check(x.c[0]), b = check(x.c[1]);
	assert(abs(a - b) <= 1);
	assert(x.h == max(a, b) + 1);
	assert(x.s == PAVL::gs(x.c[0]) + PAVL::gs(x.c[1]) + 1);
	return x.h;
}
vi vec(PAVL t) { vi r; collect(t.root, r); return r; }

PAVL build(const vi& v) {
	PAVL t;
	for (int x : v) t = t + PAVL(Node(x));
	return t;
}

void reset() { PAVL::N.clear(); PAVL::T = 0; }

// Random persistent operations on a pool of frozen versions.
void test(int ops, int maxLen) {
	reset();
	vector<PAVL> ver; vector<vi> exp; // frozen versions and their contents
	ver.emplace_back(); exp.emplace_back();
	rep(it,0,ops) {
		int i = ri(0, sz(ver) - 1), j = ri(0, sz(ver) - 1), op = ri(0, 5);
		PAVL res; vi e;
		if (op == 0) { // new tree built from scratch in this epoch
			e.resize(ri(0, 6));
			for (int& x : e) x = ri(-100, 100);
			res = build(e);
		} else if (op == 1) { // concatenate two versions (possibly the same)
			if (sz(exp[i]) + sz(exp[j]) > maxLen) continue;
			res = ver[i] + ver[j];
			e = exp[i], e.insert(e.end(), all(exp[j]));
		} else if (op == 2) { // split, keep both halves as versions
			int k = ri(0, sz(exp[i]));
			auto [l, r] = ver[i].split(k);
			assert(vec(l) == vi(exp[i].begin(), exp[i].begin() + k));
			check(l.root);
			ver.push_back(l), exp.emplace_back(exp[i].begin(), exp[i].begin() + k);
			res = r, e = vi(exp[i].begin() + k, exp[i].end());
		} else if (op == 3) { // insert a fresh node at position k
			if (sz(exp[i]) >= maxLen) continue;
			int k = ri(0, sz(exp[i])), v = ri(-100, 100);
			auto [l, r] = ver[i].split(k);
			res = l + PAVL(Node(v)) + r;
			e = exp[i], e.insert(e.begin() + k, v);
		} else if (op == 4) { // erase [a, b)
			int a = ri(0, sz(exp[i])), b = ri(a, sz(exp[i]));
			auto [l, r] = ver[i].split(a);
			auto [m, r2] = r.split(b - a);
			res = l + r2;
			e = exp[i], e.erase(e.begin() + a, e.begin() + b);
		} else { // copy [a, b) of version i into version j at position k
			int a = ri(0, sz(exp[i])), b = ri(a, sz(exp[i])), k = ri(0, sz(exp[j]));
			if (sz(exp[j]) + b - a > maxLen) continue;
			auto [l, r] = ver[i].split(a);
			auto [m, r2] = r.split(b - a);
			auto [x, y] = ver[j].split(k);
			res = x + m + y;
			e = exp[j], e.insert(e.begin() + k, exp[i].begin() + a, exp[i].begin() + b);
		}
		assert(vec(res) == e);
		check(res.root);
		ver.push_back(res), exp.push_back(e), res.step();
		if (it % 8 == 0 || it == ops - 1) rep(v,0,sz(ver)) { // old versions intact
			assert(vec(ver[v]) == exp[v]);
			check(ver[v].root);
		}
	}
}

int main(int argc, char**) {
	if (argc > 1) { // benchmark mode: ./a.out bench
		auto now = [] { return chrono::steady_clock::now(); };
		auto secs = [&](auto t0) { return chrono::duration<double>(now() - t0).count(); };
		for (int n : {100000, 500000}) {
			reset();
			auto t0 = now();
			PAVL t;
			rep(i,0,n) t = t + PAVL(Node(i)); // sorted appends
			cerr << "n=" << n << " appends: " << secs(t0) << " s, height " << PAVL::gh(t.root)
				<< ", nodes " << sz(PAVL::N) << endl;
			t0 = now();
			rep(i,0,n) { // persistent: every op is a new version
				t = t.step();
				int a = ri(0, n - 1), b = ri(a, n);
				auto [l, r] = t.split(a);
				auto [m, r2] = r.split(b - a);
				t = m + l + r2; // move [a, b) to the front
			}
			cerr << "n=" << n << " persistent cut+paste x n: " << secs(t0) << " s, height "
				<< PAVL::gh(t.root) << ", nodes " << sz(PAVL::N) << " ("
				<< sz(PAVL::N) * sizeof(Node) / 1000000 << " MB)" << endl;
			assert(PAVL::gs(t.root) == n);
		}
		reset();
		auto t0 = now();
		PAVL t = PAVL(Node(1));
		rep(i,0,29) t.step(), t = t + t; // doubling: 2^29 elements
		cerr << "doubling to size " << PAVL::gs(t.root) << ": " << secs(t0) << " s, nodes "
			<< sz(PAVL::N) << endl;
		return 0;
	}
	rep(it,0,20000) test(ri(1, 40), 10);
	rep(it,0,1000) test(200, 60);
	rep(it,0,10) test(1500, 3000);
	// edge cases: empty trees, single node, split at both ends
	reset();
	PAVL e;
	assert(vec(e + e).empty());
	{ auto [l, r] = e.split(0); assert(l.root == -1 && r.root == -1); }
	PAVL one = PAVL(Node(7));
	assert(vec(e + one) == vi{7});
	one.step();
	{ auto [l, r] = one.split(0); assert(vec(l).empty() && vec(r) == vi{7}); }
	{ auto [l, r] = one.split(1); assert(vec(l) == vi{7} && vec(r).empty()); }
	assert((vec(one + one) == vi{7, 7}));
	// sorted appends keep the tree balanced; doubling reaches huge sizes
	reset();
	PAVL t;
	rep(i,0,100000) t = t + PAVL(Node(i));
	assert(check(t.root) <= 25 && PAVL::gs(t.root) == 100000);
	rep(i,0,10) t.step(), t = t + t;
	assert(PAVL::gs(t.root) == 102400000 && PAVL::gh(t.root) <= 40);
	cout << "Tests passed!" << endl;
}
