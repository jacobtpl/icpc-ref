#include "../utilities/template.h"

#include "../../content/data-structures/BBST_Sequence_pointer.h"

// The header has no payload (Node::up / Node::down are empty and meant to be
// edited), so values live in a side array indexed by position in the pool.
// Tested: split, merge, lazy reverse, size/height bookkeeping, AVL balance.
typedef AVL::Node Node;
vector<Node> pool;
vector<int> val;
int used;

AVL make(int v) {
	assert(used < sz(pool));
	pool[used] = Node();
	val[used] = v;
	return AVL(&pool[used++]);
}

int check(Node* n, vi& out, bool rev) {
	if (!n) return 0;
	size_t before = out.size();
	int h0 = check(n->c[rev], out, rev ^ n->rev);
	out.push_back(val[n - pool.data()]);
	int h1 = check(n->c[!rev], out, rev ^ n->rev);
	assert(abs(h0 - h1) <= 1);
	assert(n->h == max(h0, h1) + 1);
	assert(n->s == (int)(out.size() - before));
	return n->h;
}
void verify(AVL t, const vi& model) {
	vi out;
	int h = check(t.root, out, 0);
	assert(out == model);
	assert(h <= 1.4405 * log2((double)sz(model) + 2) + 1e-9);
}
AVL cat(AVL a, AVL b, AVL c) { return AVL::merge(AVL::merge(a, b), c); }

void run(int K, int maxN, int ops, int checkEvery, mt19937_64& rng) {
	used = 0;
	vector<AVL> t(K);
	vector<vi> m(K);
	auto rnd = [&](int n) { return (int)(rng() % n); };
	rep(it,0,ops) {
		int a = rnd(K), n = sz(m[a]), op = rnd(100);
		int l = rnd(n + 1), r = rnd(n + 1);
		if (l > r) swap(l, r);
		if (rnd(8) == 0) l = 0;
		if (rnd(8) == 0) r = n;
		if (op < 30 && n < maxN && used < sz(pool)) { // insert at l
			int v = rnd(1000);
			auto [x, y] = t[a].split(l);
			t[a] = cat(x, make(v), y);
			m[a].insert(m[a].begin() + l, v);
		} else if (op < 45) { // erase [l, r)
			if (rnd(3)) r = min(r, l + 1 + rnd(3));
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			t[a] = AVL::merge(x, z);
			m[a].erase(m[a].begin() + l, m[a].begin() + r);
		} else if (op < 70) { // reverse [l, r)
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			if (y.root) y.root->do_rev();
			assert((y.root ? y.root->s : 0) == r - l);
			t[a] = cat(x, y, z);
			reverse(m[a].begin() + l, m[a].begin() + r);
		} else if (op < 90) { // move [l, r) of a to position p of b
			int b = rnd(K);
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			t[a] = AVL::merge(x, z);
			vi seg(m[a].begin() + l, m[a].begin() + r);
			m[a].erase(m[a].begin() + l, m[a].begin() + r);
			int p = rnd(sz(m[b]) + 1);
			auto [u, w] = t[b].split(p);
			t[b] = cat(u, y, w);
			m[b].insert(m[b].begin() + p, all(seg));
		} else { // degenerate splits and merges with empty trees
			auto [x, y] = t[a].split(0);
			assert(!x.root);
			auto [u, w] = y.split(n);
			assert(!w.root);
			t[a] = cat(AVL(), cat(x, u, w), AVL());
		}
		if (it % checkEvery == 0 || it == ops - 1)
			rep(i,0,K) verify(t[i], m[i]);
	}
}

int main() {
	mt19937_64 rng(4321);
	pool.resize(200000), val.resize(200000);
	{ // empty trees
		AVL a, b;
		assert(!AVL::merge(a, b).root);
		auto [x, y] = a.split(0);
		assert(!x.root && !y.root);
	}
	rep(it,0,20000) run(1 + (int)(rng() % 3), 1 + (int)(rng() % 8), 60, 1, rng);
	rep(it,0,2000) run(1 + (int)(rng() % 3), 40, 400, 1, rng);
	rep(it,0,20) run(2, 3000, 30000, 50, rng);
	rep(mode,0,3) { // appends / prepends / middle inserts, then repeated halving
		used = 0;
		AVL t; vi m;
		rep(i,0,100000) {
			if (mode == 0) t = AVL::merge(t, make(i)), m.push_back(i);
			else if (mode == 1) t = AVL::merge(make(i), t), m.insert(m.begin(), i);
			else {
				auto [x, y] = t.split(i / 2);
				t = cat(x, make(i), y);
				m.insert(m.begin() + i / 2, i);
			}
			if (mode == 2 && i > 3000) break;
		}
		verify(t, m);
		vector<AVL> parts;
		while (t.root && t.root->s > 1) {
			auto [x, y] = t.split(t.root->s / 2);
			parts.push_back(y), t = x;
		}
		while (!parts.empty()) t = AVL::merge(t, parts.back()), parts.pop_back();
		verify(t, m);
	}
	cout<<"Tests passed!"<<endl;
}
