#include "../utilities/template.h"

#include "../../content/various/IntervalCover.h"

// The [inclusive, inclusive] variant described in the header.
template<class T>
vi coverIncl(pair<T, T> G, vector<pair<T, T>> I) {
	vi S(sz(I)), R;
	iota(all(S), 0);
	sort(all(S), [&](int a, int b) { return I[a] < I[b]; });
	T cur = G.first;
	int at = 0;
	while (cur < G.second || R.empty()) { // (A)
		pair<T, int> mx = make_pair(cur, -1);
		while (at < sz(I) && I[S[at]].first <= cur) {
			mx = max(mx, make_pair(I[S[at]].second, S[at]));
			at++;
		}
		if (mx.second == -1) return {};
		cur = mx.first;
		R.push_back(mx.second);
	}
	return R;
}

mt19937 rng(2024);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

// incl = 0: [l, r) covering [G.first, G.second) (integer points suffice)
// incl = 1: [l, r] covering [G.first, G.second] over the reals
//           (checked on half-integer points)
template<class T>
bool covers(pair<T, T> G, const vector<pair<T, T>>& I, const vi& ids, bool incl) {
	if (!incl) {
		for (ll x = (ll)G.first; x < (ll)G.second; x++) {
			bool ok = 0;
			for (int i : ids) ok |= (ll)I[i].first <= x && x < (ll)I[i].second;
			if (!ok) return 0;
		}
		return 1;
	}
	for (ll x = 2 * (ll)G.first; x <= 2 * (ll)G.second; x++) {
		bool ok = 0;
		for (int i : ids) ok |= 2 * (ll)I[i].first <= x && x <= 2 * (ll)I[i].second;
		if (!ok) return 0;
	}
	return 1;
}

template<class T> void testType(int iters, T shift) {
	rep(it,0,iters) rep(incl,0,2) {
		int n = rnd(0, 9), C = rnd(1, 12);
		pair<T, T> G;
		int g1 = rnd(0, C), g2 = rnd(0, C);
		if (g1 > g2) swap(g1, g2);
		G = {T(g1 + shift), T(g2 + shift)};
		vector<pair<T, T>> I(n);
		for (auto& p : I) {
			int a = rnd(-1, C+1), b = rnd(-1, C+1);
			if (a > b) swap(a, b);
			p = {T(a + shift), T(b + shift)};
		}
		vi res = incl ? coverIncl(G, I) : cover(G, I);
		int best = -1;
		rep(m,0,1<<n) {
			int c = __builtin_popcount(m);
			if (best != -1 && c >= best) continue;
			vi ids;
			rep(i,0,n) if (m >> i & 1) ids.push_back(i);
			if (covers(G, I, ids, incl)) best = c;
		}
		if (!incl && g1 == g2) { assert(res.empty()); continue; }
		if (best == -1) { assert(res.empty()); continue; }
		assert(sz(res) == best);
		for (int i : res) assert(0 <= i && i < n);
		vi srt = res; sort(all(srt));
		assert(unique(all(srt)) == srt.end());
		assert(covers(G, I, res, incl));
	}
}

int main() {
	testType<int>(150000, 0);
	testType<int>(50000, -6);
	testType<ll>(30000, (ll)4e18);
	testType<ll>(30000, (ll)-4e18);
	testType<int>(30000, INT_MAX - 20);
	testType<int>(30000, INT_MIN + 20);
	testType<double>(30000, 0.0);
	// empty input
	assert(cover<int>({0, 5}, {}).empty());
	assert(cover<int>({5, 5}, {{0, 10}}).empty());
	assert(cover<int>({0, 5}, {{0, 5}}) == vi{0});
	// larger: compare size against a simple O(N^2) greedy
	rep(it,0,300) {
		int n = rnd(1, 300), C = rnd(1, 1000);
		vector<pii> I(n);
		for (auto& p : I) {
			int a = rnd(0, C), b = min(C, a + rnd(0, max(1, C / 10)));
			p = {a, b};
		}
		pii G = {rnd(0, C / 2), rnd(C / 2, C)};
		vi res = cover(G, I);
		int cur = G.first, cnt = 0; bool fail = 0;
		while (cur < G.second) {
			int nx = cur;
			for (auto& p : I) if (p.first <= cur) nx = max(nx, p.second);
			if (nx == cur) { fail = 1; break; }
			cur = nx; cnt++;
		}
		if (fail) assert(res.empty());
		else {
			assert(sz(res) == cnt);
			assert(covers(G, I, res, 0));
		}
	}
	cout << "Tests passed!" << endl;
}
