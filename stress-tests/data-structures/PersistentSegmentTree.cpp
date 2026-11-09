#include "../utilities/template.h"

#include "../../content/data-structures/PersistentSegmentTree.h"

static unsigned R = 1231231;
int ra() {
	R *= 574841;
	R += 14;
	return (int)(R >> 1);
}

ll brute(vector<ll>& v, int l, int r) {
	ll s = 0;
	rep(k,l,r) s += v[k];
	return s;
}

void check(int root, vector<ll>& v, int n, int l, int r) {
	assert(pst::query(root, 0, n, l, r).v == brute(v, l, r));
}

void run(int n, bool useBuild) {
	pst::reset();
	vector<ll> init(n);
	if (useBuild) for (ll& x : init) x = ra() % 21 - 10;
	vector<int> init32(all(init));
	int r0 = useBuild ? pst::build(init32, 0, n) : 0;
	vector<pair<int, vector<ll>>> vers = {{r0, init}};
	rep(it,0,3000) {
		pst::step();
		auto [root, v] = vers[ra() % sz(vers)];
		int ops = ra() % 4;
		rep(o,0,ops) {
			int l = ra() % (n + 1), r = ra() % (n + 1);
			if (l > r) swap(l, r);
			ll x = ra() % 21 - 10;
			root = pst::add(root, 0, n, l, r, x);
			rep(k,l,r) v[k] += x;
			int ql = ra() % (n + 1), qr = ra() % (n + 1);
			if (ql > qr) swap(ql, qr);
			check(root, v, n, ql, qr);
		}
		vers.emplace_back(root, v);
		// older versions must be intact
		rep(k,0,3) {
			auto& [oroot, ov] = vers[ra() % sz(vers)];
			int ql = ra() % (n + 1), qr = ra() % (n + 1);
			if (ql > qr) swap(ql, qr);
			check(oroot, ov, n, ql, qr);
		}
	}
	for (auto& [oroot, ov] : vers) check(oroot, ov, n, 0, n);
}

void sparse() {
	// coordinates up to 1e9 without build
	pst::reset();
	const int C = 1000000000;
	map<int, ll> v;
	vector<int> roots = {0};
	vector<map<int, ll>> vs = {v};
	rep(it,0,2000) {
		pst::step();
		int i = ra() % sz(roots);
		int p = ra() % 2 ? ra() % C : ra() % 50;
		ll x = ra() % 21 - 10;
		int root = pst::add(roots[i], 0, C, p, p + 1, x);
		auto nv = vs[i];
		nv[p] += x;
		roots.push_back(root), vs.push_back(nv);
		int j = ra() % sz(roots);
		int l = ra() % C, r = ra() % C;
		if (l > r) swap(l, r);
		ll s = 0;
		for (auto& [q, y] : vs[j]) if (l <= q && q < r) s += y;
		assert(pst::query(roots[j], 0, C, l, r).v == s);
		assert(pst::query(roots[j], 0, C, 0, C).v
			== pst::query(roots[j], 0, C, 0, 50).v
			+ pst::query(roots[j], 0, C, 50, C).v);
	}
}

int main() {
	rep(n,1,12) run(n, true), run(n, false);
	run(50, true);
	run(200, false);
	sparse();
	cout << "Tests passed!" << endl;
}
