#include "../utilities/template.h"
#define pb push_back

#include "../../content/data-structures/Treap-benq.h"

mt19937 rng(31337);
int rnd(ll a, ll b) { return (int)(a + (ll)(rng() % (unsigned long long)(b - a + 1))); }

void checkTree(pt x) { // sz and sum are up to date everywhere
	if (!x) return;
	prop(x);
	checkTree(x->c[0]), checkTree(x->c[1]);
	assert(x->sz == 1 + getsz(x->c[0]) + getsz(x->c[1]));
	assert(x->sum == x->val + getsum(x->c[0]) + getsum(x->c[1]));
}

// sequence container: inspos / delpos / splitsz / merge / flip / sum
void runSeq(int ops, int maxv) {
	pt t = nullptr;
	vi a;
	rep(it,0,ops) {
		int op = rnd(0, 5), n = sz(a);
		if (op == 0 || (op <= 2 && n < 4)) {
			int pos = rnd(0, n), v = rnd(-maxv, maxv);
			t = inspos(t, new Node(v), pos);
			a.insert(a.begin() + pos, v);
		} else if (op == 1) {
			if (!n) continue;
			int pos = rnd(0, n - 1);
			t = delpos(t, pos);
			a.erase(a.begin() + pos);
		} else {
			int l = rnd(0, n), r = rnd(0, n);
			if (l > r) swap(l, r);
			auto p1 = splitsz(t, l);
			auto p2 = splitsz(p1.second, r - l); // p2.first = [l, r)
			assert(getsz(p1.first) == l && getsz(p2.first) == r - l && getsz(p2.second) == n - r);
			assert(getsum(p2.first) == accumulate(a.begin() + l, a.begin() + r, 0LL));
			if (op <= 3 && p2.first) {
				p2.first->flip ^= 1;
				reverse(a.begin() + l, a.begin() + r);
			}
			if (op == 5) { // rotate: move [l, r) to the front
				t = merge(merge(p2.first, p1.first), p2.second);
				rotate(a.begin(), a.begin() + l, a.begin() + r);
			} else t = merge(p1.first, merge(p2.first, p2.second));
		}
		assert(getsz(t) == sz(a));
		assert(getsum(t) == accumulate(all(a), 0LL));
		if (it % 7 == 0 || it == ops - 1) {
			vi v;
			tour(t, v);
			assert(v == a);
			checkTree(t);
		}
	}
	delete t;
}

// ordered set: ins / del / split by value
void runSet(int ops, int lo, int hi) {
	pt t = nullptr;
	set<int> s;
	rep(it,0,ops) {
		int op = rnd(0, 3), v = rnd(lo, hi);
		if (op <= 1) t = ins(t, v), s.insert(v); // duplicates are not stored twice
		else if (op == 2) t = del(t, v), s.erase(v);
		else {
			auto p = split(t, v); // < v left, >= v right
			vi L, R, el, er;
			tour(p.first, L), tour(p.second, R);
			for (int x : s) (x < v ? el : er).push_back(x);
			assert(L == el && R == er);
			assert(getsz(p.first) == sz(el) && getsz(p.second) == sz(er));
			t = merge(p.first, p.second);
		}
		assert(getsz(t) == sz(s));
		assert(getsum(t) == accumulate(all(s), 0LL));
		if (it % 5 == 0 || it == ops - 1) {
			vi v;
			tour(t, v);
			assert(v == vi(all(s)));
			checkTree(t);
		}
	}
	delete t;
}

int main() {
	srand(5);
	{ // empty / single element
		assert(merge(nullptr, nullptr) == nullptr);
		assert(split(nullptr, 3).first == nullptr && splitsz(nullptr, 0).second == nullptr);
		assert(del(nullptr, 1) == nullptr && getsz(nullptr) == 0 && getsum(nullptr) == 0);
		vi v; tour(nullptr, v); assert(v.empty());
		pt t = ins(nullptr, 5);
		auto p = splitsz(t, 1);
		assert(p.first == t && !p.second);
		p = splitsz(t, 0);
		assert(!p.first && p.second == t);
		t->flip ^= 1;
		tour(t, v); assert(v == vi{5});
		delete t;
	}
	rep(it,0,20000) runSeq(rnd(1, 40), it % 2 ? 3 : 1000000000);
	rep(it,0,200) runSeq(1500, 1000000000);
	rep(it,0,20000) runSet(rnd(1, 40), -5, 5);
	rep(it,0,10000) runSet(rnd(1, 40), INT_MIN, INT_MIN + 6);
	rep(it,0,10000) runSet(rnd(1, 40), INT_MAX - 7, INT_MAX - 1); // ins/del use v+1
	rep(it,0,200) runSet(1500, -300, 300);
	rep(it,0,200) runSet(1500, INT_MIN, INT_MAX - 1);
	{ // sorted insertion order (worst case for an unbalanced BST), sums near ll range
		pt t = nullptr;
		int n = 200000;
		rep(i,0,n) t = inspos(t, new Node(INT_MAX - i), i);
		assert(getsz(t) == n && getsum(t) == (ll)n * INT_MAX - (ll)n * (n - 1) / 2);
		t->flip ^= 1;
		vi v; tour(t, v);
		rep(i,0,n) assert(v[i] == INT_MAX - (n - 1 - i));
		delete t;
	}
	cout<<"Tests passed!"<<endl;
}
