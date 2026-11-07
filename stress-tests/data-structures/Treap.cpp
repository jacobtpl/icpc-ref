#include "../utilities/template.h"

#include "../../content/data-structures/Treap.h"

vi tour(pt t) { vi v; each(t, [&](int x) { v.push_back(x); }); return v; }

int main() {
	srand(7);
	// sequential container: inspos/delpos/rev/range sum
	rep(it,0,300) {
		pt t = nullptr; vi ref;
		rep(op,0,300) {
			int k = rand() % 4, n = sz(ref);
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
					t = rev(t, l, r);
					reverse(ref.begin() + l, ref.begin() + r);
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
		pt t = nullptr; set<int> ref;
		rep(op,0,300) {
			int v = rand() % 50;
			if (rand() % 2) t = ins(t, v), ref.insert(v);
			else t = del(t, v), ref.erase(v);
			assert(tour(t) == vi(all(ref)));
			int q = rand() % 50;
			auto p = split(t, q);
			assert(getsz(p.first) == (int)distance(ref.begin(), ref.lower_bound(q)));
			t = merge(p.first, p.second);
		}
	}
	cout<<"Tests passed!"<<endl;
}
