#include "../utilities/template.h"

#include "../../content/data-structures/Treap.h"

typedef vector<ll> vl;
vl tour(pt t) { vl v; each(t, [&](ll x) { v.push_back(x); }); return v; }

int main() {
	srand(7);
	// sequential container: inspos/delpos/rev/range sum
	rep(it,0,300) {
		vl ref(rand() % 20);
		for (ll& x : ref) x = rand() % 100;
		pt t = build(ref);
		assert(tour(t) == ref);
		rep(op,0,300) {
			int k = rand() % 6, n = sz(ref);
			if (k == 0 || ref.empty()) {
				int pos = rand() % (n + 1), v = rand() % 100;
				t = inspos(t, new Node(v), pos);
				ref.insert(ref.begin() + pos, v);
			} else if (k == 1) {
				int pos = rand() % n;
				t = delpos(t, pos);
				ref.erase(ref.begin() + pos);
			} else {
				int l = rand() % n, r = rand() % n;
				if (l > r) swap(l, r);
				r++;
				if (k == 2) {
					t = onRange(t, l, r, [&](pt x) { x->flip ^= 1; });
					reverse(ref.begin() + l, ref.begin() + r);
				} else if (k == 3) {
					ll v = rand() % 21 - 10;
					t = onRange(t, l, r, [&](pt x) { upd(x, v); });
					rep(i,l,r) ref[i] += v;
				} else if (k == 4) {
					ll mn = LLONG_MAX;
					t = onRange(t, l, r, [&](pt x) { mn = x->mn; });
					assert(mn == *min_element(ref.begin() + l, ref.begin() + r));
				} else {
					auto a = splitsz(t, l), b = splitsz(a.second, r - l);
					ll s = 0;
					rep(i,l,r) s += ref[i];
					assert(getsum(b.first) == s);
					assert(getsz(b.first) == r - l);
					t = merge(a.first, merge(b.first, b.second));
				}
			}
			assert(tour(t) == ref);
		}
	}
	// ordered set: ins/del/split by key
	rep(it,0,300) {
		pt t = nullptr; set<ll> ref;
		rep(op,0,300) {
			ll v = rand() % 50;
			if (rand() % 2) t = ins(t, v), ref.insert(v);
			else t = del(t, v), ref.erase(v);
			assert(tour(t) == vl(all(ref)));
			int q = rand() % 50;
			auto p = split(t, q);
			assert(getsz(p.first) == (int)distance(ref.begin(), ref.lower_bound(q)));
			t = merge(p.first, p.second);
		}
	}
	cout<<"Tests passed!"<<endl;
}
