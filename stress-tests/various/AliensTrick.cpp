#include "../utilities/template.h"

#include "../../content/various/AliensTrick.h"

typedef __int128 LL;

mt19937_64 rng(1);
ll rnd(ll lo, ll hi) { return (ll)(rng() % (unsigned long long)(hi - lo + 1)) + lo; }

int calls;
// mode 0: smallest optimal count, 1: largest, 2: random optimal count
template<class T>
auto solver(const vector<T>& c, int mode) {
	return [&c, mode](ll lam) {
		calls++;
		LL best = LL(c[0]);
		rep(j,1,sz(c)) best = min(best, LL(c[j]) + LL(lam) * j);
		vi opt;
		rep(j,0,sz(c)) if (LL(c[j]) + LL(lam) * j == best) opt.push_back(j);
		int cnt = mode == 0 ? opt[0] : mode == 1 ? opt.back()
			: opt[rnd(0, sz(opt) - 1)];
		assert(best >= LLONG_MIN && best <= LLONG_MAX);
		return pair<ll, int>((ll)best, cnt);
	};
}

int maxCalls(ll lo, ll hi) {
	int r = 0; // floor(log2(len)) + 1
	for (LL len = LL(hi) - lo + 1; len >= 1; len /= 2) r++;
	return r;
}

// convex table c[0..n] with integer slopes in [-S, S], values near base
vector<ll> convexTable(int n, ll S, ll base) {
	vector<ll> d(n);
	for (ll& x : d) x = rnd(-S, S);
	sort(all(d));
	vector<ll> c(n + 1, base);
	rep(i,0,n) c[i + 1] = c[i] + d[i];
	return c;
}

void testTable(const vector<ll>& c, ll S, bool wide, bool checkCalls) {
	int n = sz(c) - 1;
	rep(k,0,n + 1) rep(mode,0,3) {
		// exact: lambda in [-(c[k+1]-c[k]), -(c[k]-c[k-1])] makes k optimal
		ll lo = -S, hi = S;
		if (!wide) {
			lo = k < n ? -(c[k + 1] - c[k]) : -S;
			hi = k > 0 ? -(c[k] - c[k - 1]) : S;
			if (lo > hi) swap(lo, hi);
		}
		calls = 0;
		ll got = aliens(lo, hi, k, solver(c, mode));
		if (got != c[k]) {
			cout << "FAIL table n=" << n << " k=" << k << " mode=" << mode
				 << " lo=" << lo << " hi=" << hi << " got " << got
				 << " want " << c[k] << endl;
			for (ll x : c) cout << x << ' ';
			cout << endl;
			abort();
		}
		if (checkCalls) assert(calls <= maxCalls(lo, hi));
	}
}

void testTables() {
	rep(it,0,100000) {
		int n = (int)rnd(0, 8);
		ll S = rnd(0, 1) ? rnd(0, 5) : rnd(0, (ll)1e12);
		ll base = rnd(-(ll)1e12, (ll)1e12);
		vector<ll> c = convexTable(n, S, base);
		testTable(c, S, it & 1, true);
	}
	// linear c and all-equal slopes, duplicates
	rep(it,0,2000) {
		int n = (int)rnd(1, 8);
		ll d = rnd(-10, 10);
		vector<ll> c(n + 1);
		rep(i,0,n + 1) c[i] = 7 + d * i;
		testTable(c, 10, it & 1, true);
	}
	// near overflow: slopes up to 1e18, values up to ~4e18/8 so lam*k fits
	rep(it,0,3000) {
		int n = (int)rnd(0, 4);
		ll S = rnd((ll)1e17, (ll)1e18);
		vector<ll> c = convexTable(n, S, rnd(-(ll)1e17, (ll)1e17));
		testTable(c, S, it & 1, true);
	}
}

// pick exactly k pairwise non-adjacent elements with minimum sum
void testNonAdjacent() {
	rep(it,0,4000) {
		int n = (int)rnd(0, 10);
		ll V = rnd(0, 1) ? 3 : (ll)1e15;
		vector<ll> a(n);
		for (ll& x : a) x = rnd(-V, V);
		if (rnd(0, 3) == 0) fill(all(a), rnd(-V, V));
		int K = (n + 1) / 2;
		// dp[i][j] = min sum picking j from first i, non-adjacent
		const ll INF = LLONG_MAX / 4;
		vector<vector<ll>> dp(n + 1, vector<ll>(K + 1, INF));
		dp[0][0] = 0;
		rep(i,0,n) rep(j,0,K + 1) if (dp[i][j] < INF) {
			dp[i + 1][j] = min(dp[i + 1][j], dp[i][j]);
			if (j < K) {
				ll prev = i ? dp[i - 1][j] : (j == 0 ? 0 : INF);
				if (i == 0 && j > 0) prev = INF;
				if (prev < INF) dp[i + 1][j + 1] = min(dp[i + 1][j + 1], prev + a[i]);
			}
		}
		ll C = 0; // documented bound: max |c(j) - c(j-1)|
		rep(j,1,K + 1) C = max(C, abs(dp[n][j] - dp[n][j - 1]));
		rep(k,0,K + 1) rep(mode,0,3) {
			auto f = [&](ll lam) {
				calls++;
				// O(n) dp: (value, count) with tie-breaking by mode
				typedef pair<ll, int> P;
				auto better = [&](P x, P y) {
					if (x.first != y.first) return x.first < y.first;
					if (mode == 0) return x.second < y.second;
					if (mode == 1) return x.second > y.second;
					return (bool)(rng() & 1);
				};
				P take(0, 0), skip(0, 0), prev(0, 0); // dp over prefix
				rep(i,0,n) {
					P t(prev.first + a[i] + lam, prev.second + 1);
					P s = better(take, skip) ? take : skip;
					prev = s;
					take = t; skip = s;
				}
				return better(take, skip) ? take : skip;
			};
			ll got = aliens(-C, C, k, f);
			if (got != dp[n][k]) {
				cout << "FAIL nonadj n=" << n << " k=" << k << " mode=" << mode
					 << " got " << got << " want " << dp[n][k] << endl;
				for (ll x : a) cout << x << ' ';
				cout << endl;
				abort();
			}
		}
	}
}

// split into exactly k segments minimising sum of squared segment sums
void testSegments() {
	rep(it,0,2000) {
		int n = (int)rnd(1, 8);
		vector<ll> a(n), pre(n + 1, 0);
		for (ll& x : a) x = rnd(0, 6); // non-negative: c(k) convex (Monge cost)
		rep(i,0,n) pre[i + 1] = pre[i] + a[i];
		auto cost = [&](int l, int r) { ll s = pre[r] - pre[l]; return s * s; };
		const ll INF = LLONG_MAX / 4;
		vector<vector<ll>> dp(n + 1, vector<ll>(n + 1, INF));
		dp[0][0] = 0;
		rep(i,1,n + 1) rep(j,1,i + 1) rep(p,j - 1,i)
			if (dp[p][j - 1] < INF) dp[i][j] = min(dp[i][j], dp[p][j - 1] + cost(p, i));
		ll C = 0;
		rep(j,2,n + 1) C = max(C, abs(dp[n][j] - dp[n][j - 1]));
		rep(k,1,n + 1) rep(mode,0,3) {
			auto f = [&](ll lam) {
				calls++;
				typedef pair<ll, int> P;
				vector<P> d(n + 1, P(INF, 0));
				d[0] = P(0, 0);
				rep(i,1,n + 1) rep(p,0,i) {
					P cand(d[p].first + cost(p, i) + lam, d[p].second + 1);
					bool b = cand.first < d[i].first || (cand.first == d[i].first &&
						(mode == 0 ? cand.second < d[i].second
						: mode == 1 ? cand.second > d[i].second : (bool)(rng() & 1)));
					if (b) d[i] = cand;
				}
				return d[n];
			};
			ll got = aliens(-C, C, k, f);
			if (got != dp[n][k]) {
				cout << "FAIL segments n=" << n << " k=" << k << " mode=" << mode
					 << " got " << got << " want " << dp[n][k] << endl;
				abort();
			}
		}
	}
}

int main() {
	testTables();
	testNonAdjacent();
	testSegments();
	cout << "Tests passed!" << endl;
}
