#include "../utilities/template.h"

#include "../../content/data-structures/PBBST.h"

void tour(int n, vi& v) {
	if (n == -1) return;
	tour(PAVL::N[n].c[0], v);
	v.push_back(PAVL::N[n].val);
	tour(PAVL::N[n].c[1], v);
}

int main() {
	srand(5);
	rep(it,0,200) {
		PAVL::N.clear(); PAVL::T = 0;
		vector<PAVL> vers; vector<vi> refs;
		PAVL cur; vi ref;
		rep(op,0,200) {
			cur = cur.step();
			int k = rand() % 4, n = sz(ref);
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
			} else if (k == 2) { // move [l, r) to position kk
				int l = rand() % n, r = rand() % n;
				if (l > r) swap(l, r);
				r++;
				int kk = rand() % (n - (r - l) + 1);
				auto [a, b] = cur.split(l);
				auto [m, c] = b.split(r - l);
				PAVL rest = a + c;
				auto [x, y] = rest.split(kk);
				cur = x + m + y;
				vi seg(ref.begin() + l, ref.begin() + r);
				ref.erase(ref.begin() + l, ref.begin() + r);
				ref.insert(ref.begin() + kk, all(seg));
			} else if (!vers.empty()) { // go back to an old version
				int j = rand() % sz(vers);
				cur = vers[j]; ref = refs[j];
			}
			vers.push_back(cur); refs.push_back(ref);
			rep(j,0,sz(vers)) {
				vi v; tour(vers[j].root, v);
				assert(v == refs[j]);
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
