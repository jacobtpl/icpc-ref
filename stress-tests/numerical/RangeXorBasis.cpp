#include "../utilities/template.h"

#include "../../content/numerical/RangeXorBasis.h"

int main() {
	rep(it,0,30000) {
		int n = rand() % 12 + 1, bits = rand() % 6 + 1;
		if (it % 10 == 0) bits = 30;
		vector<unsigned> a(n);
		for (auto &x : a) x = bits == 30 ? (unsigned)rand() & ((1u << 30) - 1) : rand() % (1 << bits);
		if (bits != 30 && rand() % 3 == 0) for (auto &x : a) if (rand() % 2) x = 0;
		RangeXorBasis rb;
		for (auto x : a) rb.push(x);
		rep(l,0,n) rep(r,l,n) {
			if (bits <= 6) {
				vector<char> span(1 << bits); span[0] = 1; int cnt = 1;
				rep(i,l,r+1) rep(v,0,1<<bits) if (span[v] && !span[v ^ a[i]]) {
					// closure: add a[i] to every current element
					vector<char> ns = span;
					rep(u,0,1<<bits) if (span[u]) ns[u ^ a[i]] = 1;
					span = ns; break;
				}
				cnt = (int)count(all(span), 1);
				assert(cnt == 1 << rb.rank(l, r));
				unsigned x = rand() % (1 << bits), mx = 0, mn = ~0u;
				rep(v,0,1<<bits) if (span[v]) mx = max(mx, x ^ v), mn = min(mn, x ^ v);
				assert(rb.maxXor(l, r, x) == mx);
				assert(rb.minXor(l, r, x) == mn);
				assert(rb.has(l, r, x) == (bool)span[x]);
			} else {
				unsigned mx = 0;
				rep(m,0,1<<(r-l+1)) { unsigned v = 0; rep(i,l,r+1) if (m >> (i-l) & 1) v ^= a[i]; mx = max(mx, v); }
				assert(rb.maxXor(l, r) == mx);
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
