#include "../utilities/template.h"

#include "../../content/strings/ParamAhoCorasick.h"

bool iso(const vi& a, const vi& b) {
	map<int,int> f, g;
	rep(i,0,sz(a)) {
		if (f.count(a[i]) && f[a[i]] != b[i]) return 0;
		if (g.count(b[i]) && g[b[i]] != a[i]) return 0;
		f[a[i]] = b[i]; g[b[i]] = a[i];
	}
	return 1;
}

int main() {
	rep(it,0,30000) {
		int al = rand() % 4 + 1, k = rand() % 6 + 1, n = rand() % 30;
		if (it % 100 == 0) al = 30, n = 200;
		vector<vi> pat(k); vector<ll> w(k);
		ParamAC ac;
		rep(i,0,k) {
			pat[i].resize(rand() % 8 + 1);
			for (int &x : pat[i]) x = rand() % al;
			w[i] = rand() % 100 + 1;
			ac.add(pat[i], w[i]);
		}
		ac.build();
		vi s(n);
		for (int &x : s) x = rand() % al;
		vector<ll> got(n);
		ac.run(s, [&](int i, int v) { got[i] = ac.t[v].val; });
		rep(i,0,n) {
			ll want = 0;
			rep(j,0,k) if (sz(pat[j]) <= i + 1 && iso(pat[j],
				vi(s.begin() + i + 1 - sz(pat[j]), s.begin() + i + 1)))
				want += w[j];
			assert(got[i] == want);
		}
	}
	// strings, UCup S5 Stage 8 E sample: 3 windows of "AABBCC" match "aa"
	ParamAC ac; ac.add(string("aa")); ac.build();
	ll c = 0;
	ac.run(string("AABBCC"), [&](int, int v) { c += ac.t[v].val; });
	assert(c == 3);
	cout << "Tests passed!" << endl;
}
