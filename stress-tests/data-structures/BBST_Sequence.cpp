#include "../utilities/template.h"

#include "../../content/data-structures/BBST_Sequence.h"

// Payload is kept in side arrays indexed by node id, since the header leaves
// Node::up / Node::down to the user: value, subtree sum, lazy add, lazy reverse.
vector<ll> val, sum, ad;
vector<char> rv;
Node* base; // == AVL::N.data()

void applyAdd(int i, ll x) { val[i] += x, sum[i] += x * base[i].s, ad[i] += x; }
void applyRev(int i) { rv[i] ^= 1, swap(base[i].c[0], base[i].c[1]); }
void Node::up() {
	base = AVL::N.data();
	int i = (int)(this - base);
	sum[i] = val[i];
	for (int ch : c) if (ch != -1) sum[i] += sum[ch];
}
void Node::down() {
	base = AVL::N.data();
	int i = (int)(this - base);
	for (int ch : c) if (ch != -1) {
		if (ad[i]) applyAdd(ch, ad[i]);
		if (rv[i]) applyRev(ch);
	}
	ad[i] = 0, rv[i] = 0;
}

int live = 0;
AVL make(ll v) {
	AVL t = AVL::make_avl(Node());
	if (t.root >= sz(val)) val.resize(t.root + 1), sum.resize(t.root + 1), ad.resize(t.root + 1), rv.resize(t.root + 1);
	val[t.root] = sum[t.root] = v, ad[t.root] = 0, rv[t.root] = 0;
	base = &t.get_root() - t.root;
	live++;
	return t;
}

// Checks size/height/sum bookkeeping and the AVL balance invariant, and
// appends the in-order contents (lazies are resolved on the way, not pushed).
int check(int n, vector<ll>& out, ll add, bool rev) {
	if (n == -1) return 0;
	Node& x = base[n];
	size_t before = out.size();
	ll cadd = add + ad[n]; bool crev = rev ^ rv[n];
	int h0 = check(x.c[rev], out, cadd, crev);
	size_t mid = out.size();
	out.push_back(val[n] + add);
	int h1 = check(x.c[!rev], out, cadd, crev);
	assert(abs(h0 - h1) <= 1);
	assert(x.h == max(h0, h1) + 1);
	assert(x.s == (int)(out.size() - before));
	ll s = 0;
	rep(i,(int)before,sz(out)) s += out[i];
	assert(sum[n] + add * x.s == s);
	(void)mid;
	return x.h;
}
void verify(AVL& t, const vector<ll>& model) {
	vector<ll> out;
	int h = check(t.root, out, 0, 0);
	assert(out == model);
	// AVL height bound: h <= 1.4405 log2(n + 2)
	assert(h <= 1.4405 * log2((double)sz(model) + 2) + 1e-9);
}

void run(int K, int maxN, int ops, int checkEvery, ll V, mt19937_64& rng) {
	vector<AVL> t(K);
	vector<vector<ll>> m(K);
	auto rnd = [&](int n) { return (int)(rng() % n); };
	auto rv_ = [&]() { return (ll)(rng() % (2 * V + 1)) - V; };
	rep(it,0,ops) {
		int a = rnd(K), n = sz(m[a]), op = rnd(100);
		int l = rnd(n + 1), r = rnd(n + 1);
		if (l > r) swap(l, r);
		if (rnd(8) == 0) l = 0;
		if (rnd(8) == 0) r = n;
		if (op < 30 && n < maxN) { // insert at l
			ll v = rv_();
			auto [x, y] = t[a].split(l);
			t[a] = x + make(v) + y;
			m[a].insert(m[a].begin() + l, v);
		} else if (op < 45) { // erase [l, r)
			if (rnd(3)) r = min(r, l + 1 + rnd(3));
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			live -= r - l;
			y.clear();
			assert(y.root == -1);
			t[a] = x + z;
			m[a].erase(m[a].begin() + l, m[a].begin() + r);
		} else if (op < 60) { // reverse [l, r)
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			if (y.root != -1) applyRev(y.root);
			t[a] = x + y + z;
			reverse(m[a].begin() + l, m[a].begin() + r);
		} else if (op < 72) { // add on [l, r)
			ll v = rv_();
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			if (y.root != -1) applyAdd(y.root, v);
			t[a] = x + y + z;
			rep(i,l,r) m[a][i] += v;
		} else if (op < 85) { // sum of [l, r)
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			ll s = 0;
			rep(i,l,r) s += m[a][i];
			if (y.root == -1) assert(l == r);
			else assert(sum[y.root] == s && y.get_root().s == r - l);
			t[a] = x + y + z;
		} else if (op < 95) { // move [l, r) of a to position p of b
			int b = rnd(K);
			auto [x, yz] = t[a].split(l);
			auto [y, z] = yz.split(r - l);
			t[a] = x + z;
			vector<ll> seg(m[a].begin() + l, m[a].begin() + r);
			m[a].erase(m[a].begin() + l, m[a].begin() + r);
			int p = rnd(sz(m[b]) + 1);
			auto [u, w] = t[b].split(p);
			u += y; u += w;
			t[b] = u;
			m[b].insert(m[b].begin() + p, all(seg));
		} else if (op < 97) { // clear everything
			live -= n;
			t[a].clear();
			m[a].clear();
		} else { // degenerate splits and merges with empty trees
			auto [x, y] = t[a].split(0);
			assert(x.root == -1);
			auto [u, w] = y.split(n);
			assert(w.root == -1);
			t[a] = AVL() + (x + u) + w + AVL();
		}
		if (it % checkEvery == 0 || it == ops - 1)
			rep(i,0,K) verify(t[i], m[i]);
	}
	rep(i,0,K) live -= sz(m[i]), t[i].clear();
	assert(live == 0);
}

int main() {
	mt19937_64 rng(1234);
	{ // empty trees
		AVL a, b;
		a += b;
		assert(a.root == -1 && (a + b).root == -1);
		auto [x, y] = a.split(0);
		assert(x.root == -1 && y.root == -1);
		a.clear();
		AVL c = make(5);
		assert(c.get_root().s == 1 && sum[c.root] == 5);
		c.clear(); live--;
	}
	rep(it,0,20000) run(1 + (int)(rng() % 3), 1 + (int)(rng() % 8), 60, 1, 5, rng);
	rep(it,0,2000) run(1 + (int)(rng() % 3), 40, 400, 1, 1000000000, rng);
	rep(it,0,20) run(2, 3000, 30000, 50, 1000000, rng);
	{ // sequential appends / prepends / middle inserts, then repeated halving
		rep(mode,0,3) {
			AVL t; vector<ll> m;
			rep(i,0,100000) {
				if (mode == 0) t += make(i), m.push_back(i);
				else if (mode == 1) t = make(i) + t, m.insert(m.begin(), i);
				else {
					auto [x, y] = t.split(i / 2);
					t = x + make(i) + y;
					m.insert(m.begin() + i / 2, i);
				}
				if (mode == 2 && i > 3000) break;
			}
			verify(t, m);
			vector<AVL> parts;
			while (t.root != -1 && t.get_root().s > 1) {
				auto [x, y] = t.split(t.get_root().s / 2);
				parts.push_back(y), t = x;
			}
			while (!parts.empty()) t += parts.back(), parts.pop_back();
			verify(t, m);
			live -= sz(m), t.clear();
		}
		assert(live == 0);
	}
	cout<<"Tests passed!"<<endl;
}
