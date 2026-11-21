#include "../utilities/template.h"
#include "../utilities/random.h"
#define pb push_back

#include "../../content/graph/Rerooting.h"

// Naive: value that u sends to its neighbour g[u][i] (i = -1: answer for u as root).
template<class T, class M, class F>
T naiveSend(const vector<vi>& g, int u, int i, T e, M& mg, F& f) {
	T x = e;
	rep(j,0,sz(g[u])) if (j != i) {
		int w = g[u][j];
		int back = int(find(all(g[w]), u) - g[w].begin());
		x = mg(x, naiveSend(g, w, back, e, mg, f));
	}
	return f(x, u, i);
}

template<class T, class M, class F>
vector<T> naive(const vector<vi>& g, T e, M mg, F f) {
	vector<T> res(sz(g));
	rep(r,0,sz(g)) res[r] = naiveSend(g, r, -1, e, mg, f);
	return res;
}

vector<vi> genForest(int n, bool forest) {
	int shape = rand() % 5;
	vector<pii> ed;
	rep(v,1,n) {
		int p;
		if (shape == 0) p = rand() % v;
		else if (shape == 1) p = v - 1;
		else if (shape == 2) p = 0;
		else if (shape == 3) p = max(0, v - 1 - rand() % 3);
		else p = (v - 1) / 2;
		if (forest && rand() % 3 == 0) continue;
		ed.pb({p, v});
	}
	vi perm(n);
	iota(all(perm), 0);
	shuffle_vec(perm);
	vector<vi> g(n);
	for (auto [a, b] : ed) {
		g[perm[a]].pb(perm[b]);
		g[perm[b]].pb(perm[a]);
	}
	for (auto& adj : g) shuffle_vec(adj);
	return g;
}

typedef array<array<int, 2>, 2> Mat;
Mat matMul(Mat a, Mat b) {
	Mat c{};
	rep(i,0,2) rep(j,0,2) rep(k,0,2) c[i][j] = (c[i][j] + a[i][k] * b[k][j]) % 7;
	return c;
}

void testCase(const vector<vi>& g) {
	int n = sz(g);

	// 1. Exact merge order and (v, i) arguments via string concatenation.
	auto cat = [](string a, string b) { return a + b; };
	auto rec = [](string x, int v, int i) {
		return "(" + x + ")" + to_string(v) + "," + to_string(i) + ";";
	};
	assert(reroot(g, string(), cat, rec) == naive(g, string(), cat, rec));

	// 2. Longest weighted path, negative weights and values near +-1e18.
	vector<vector<ll>> w(n);
	bool big = randBool(), nonneg = randBool();
	map<pii, ll> ew;
	rep(v,0,n) for (int u : g[v]) {
		pii key = minmax(u, v);
		if (!ew.count(key)) {
			ll x = big ? (ll)randRange(uint64_t(1e18 / max(n, 1))) : (ll)randRange(20);
			if (!nonneg && randBool()) x = -x;
			ew[key] = x;
		}
		w[v].pb(ew[key]);
	}
	auto mx = [](ll a, ll b) { return max(a, b); };
	auto far = [&](ll x, int v, int i) { return i < 0 ? x : x + w[v][i]; };
	vector<ll> F = reroot(g, 0LL, mx, far);
	assert(F == naive(g, 0LL, mx, far));
	if (nonneg) rep(r,0,n) { // independent oracle: DFS from r
		vector<ll> d(n, -1);
		vi st = {r};
		d[r] = 0;
		while (!st.empty()) {
			int v = st.back(); st.pop_back();
			rep(i,0,sz(g[v])) if (d[g[v][i]] < 0)
				d[g[v][i]] = d[v] + w[v][i], st.pb(g[v][i]);
		}
		assert(F[r] == *max_element(all(d)));
	}

	// 3. (subtree size, sum of distances) vs BFS.
	typedef pair<ll, ll> P;
	auto add = [](P a, P b) { return P(a.first + b.first, a.second + b.second); };
	auto sub = [](P x, int, int i) { return i < 0 ? x : P(x.first + 1, x.second + x.first + 1); };
	vector<P> S = reroot(g, P(0, 0), add, sub);
	rep(r,0,n) {
		vector<ll> d(n, -1);
		vi q = {r};
		d[r] = 0;
		rep(j,0,sz(q)) for (int u : g[q[j]]) if (d[u] < 0) d[u] = d[q[j]] + 1, q.pb(u);
		ll cnt = 0, sum = 0;
		for (ll x : d) if (x >= 0) cnt++, sum += x;
		assert(S[r] == P(cnt - 1, sum));
	}

	// 4. Non-commutative, non-invertible monoid: 2x2 matrices mod 7.
	vector<vector<Mat>> em(n);
	rep(v,0,n) rep(i,0,sz(g[v])) {
		Mat m;
		rep(a,0,2) rep(b,0,2) m[a][b] = rand() % 4 == 0 ? 0 : rand() % 7;
		em[v].pb(m);
	}
	Mat I{{{1, 0}, {0, 1}}};
	auto mf = [&](Mat x, int v, int i) { return i < 0 ? x : matMul(em[v][i], x); };
	assert(reroot(g, I, matMul, mf) == naive(g, I, matMul, mf));
}

int main() {
	testCase(vector<vi>());
	testCase(vector<vi>(1));
	testCase(vector<vi>(5)); // forest of isolated vertices
	rep(it,0,20000) testCase(genForest(rand() % 10, rand() % 4 == 0));
	rep(it,0,200) testCase(genForest(10 + rand() % 30, randBool()));

	// Large path: non-recursive, O(N); answers checked against closed form.
	int n = 1000000;
	vector<vi> g(n);
	rep(v,1,n) g[v - 1].pb(v), g[v].pb(v - 1);
	typedef pair<ll, ll> P;
	vector<P> S = reroot(g, P(0, 0), [](P a, P b) { return P(a.first + b.first, a.second + b.second); },
		[](P x, int, int i) { return i < 0 ? x : P(x.first + 1, x.second + x.first + 1); });
	rep(r,0,n) {
		ll a = r, b = n - 1 - r;
		assert(S[r] == P(n - 1, a * (a + 1) / 2 + b * (b + 1) / 2));
	}
	cout << "Tests passed!" << endl;
}
