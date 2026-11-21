#include "../utilities/template.h"

#include "../../content/data-structures/PBBST.h"

void tour(int n, vector<ll>& v, ll add) { // read-only, no copies
	if (n == -1) return;
	auto& x = PAVL::N[n];
	tour(x.c[0], v, add + x.add);
	v.push_back(x.val + add);
	tour(x.c[1], v, add + x.add);
}

int main() {
	srand(5);
	rep(it,0,300) {
		PAVL::reset();
		vector<PAVL> vers; vector<vector<ll>> refs;
		PAVL cur; vector<ll> ref;
		rep(op,0,150) {
			rep(sub,0,1 + rand() % 3) { // several ops per version
				int k = rand() % 7, n = sz(ref);
				if (k == 0 || ref.empty()) {
					int pos = rand() % (n + 1), v = rand() % 100;
					auto [l, r] = cur.split(pos);
					cur = l + PAVL(Node(v)) + r;
					ref.insert(ref.begin() + pos, v);
				} else if (k == 1) {
					int pos = rand() % n;
					auto [l, r] = cur.split(pos);
					auto [m, rr] = r.split(1);
					cur = l + rr;
					ref.erase(ref.begin() + pos);
				} else if (k == 6) { // go back to an old version
					if (vers.empty()) continue;
					int j = rand() % sz(vers);
					cur = vers[j]; ref = refs[j];
				} else {
					int l = rand() % n, r = rand() % n;
					if (l > r) swap(l, r);
					r++;
					auto [a, b] = cur.split(l);
					auto [m, c] = b.split(r - l);
					if (k == 2) { // move [l, r) to position kk
						int kk = rand() % (n - (r - l) + 1);
						PAVL rest = a + c;
						auto [x, y] = rest.split(kk);
						cur = x + m + y;
						vector<ll> seg(ref.begin() + l, ref.begin() + r);
						ref.erase(ref.begin() + l, ref.begin() + r);
						ref.insert(ref.begin() + kk, all(seg));
						continue;
					}
					if (k == 3) {
						int v = rand() % 100 - 50;
						m.get_root().apply(v);
						rep(i,l,r) ref[i] += v;
					} else {
						assert(m.get_root().sum == accumulate(ref.begin() + l, ref.begin() + r, 0LL));
						assert(m.get_root().s == r - l);
					}
					cur = a + m + c;
				}
			}
			if (rand() % 4 == 0 && cur.root != -1) { // whole-tree add, no split
				cur.get_root().apply(7);
				for (ll& x : ref) x += 7;
			}
			vers.push_back(cur); refs.push_back(ref);
			PAVL::step();
			rep(j,0,sz(vers)) {
				vector<ll> v; tour(vers[j].root, v, 0);
				assert(v == refs[j]);
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
