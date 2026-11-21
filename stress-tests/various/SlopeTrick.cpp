#include "../utilities/template.h"

#include "../../content/various/SlopeTrick.h"

// Dense oracle: f on [-B, B], linear outside (all breakpoints stay inside).
const int B = 60, V = 10, S = 4;
struct Naive {
	vector<ll> f = vector<ll>(2 * B + 1, 0);
	ll& at(ll x) { return f[x + B]; }
	ll eval(ll x) {
		if (x < -B) return at(-B) + (at(-B) - at(-B + 1)) * (-B - x);
		if (x > B) return at(B) + (at(B) - at(B - 1)) * (x - B);
		return at(x);
	}
	void add(ll c) { for (ll& v : f) v += c; }
	void addR(ll a) { rep(x,-B,B+1) at(x) += max(x - a, 0LL); }
	void addL(ll a) { rep(x,-B,B+1) at(x) += max(a - x, 0LL); }
	void addAbs(ll a) { addR(a), addL(a); }
	void prefMin() { rep(x,-B+1,B+1) at(x) = min(at(x), at(x - 1)); }
	void sufMin() { for (ll x = B - 1; x >= -B; x--) at(x) = min(at(x), at(x + 1)); }
	void shift(ll a, ll b) {
		vector<ll> g(2 * B + 1);
		rep(x,-B,B+1) {
			ll best = LLONG_MAX;
			for (ll y = x - b; y <= x - a; y++) best = min(best, eval(y));
			g[x + B] = best;
		}
		f = g;
	}
	ll mn() { return *min_element(all(f)); }
	ll topL() { return at(-B) == mn() ? -SlopeTrick::inf : find(all(f), mn()) - f.begin() - B; }
	ll topR() { return at(B) == mn() ? SlopeTrick::inf : B - (find(f.rbegin(), f.rend(), mn()) - f.rbegin()); }
};

ll evalST(SlopeTrick& st, ll x) {
	ll r = st.mn;
	auto L = st.L; while (sz(L)) r += max(0LL, L.top() + st.aL - x), L.pop();
	auto R = st.R; while (sz(R)) r += max(0LL, x - (R.top() + st.aR)), R.pop();
	return r;
}

void check(SlopeTrick& st, Naive& nv, SlopeTrick& far, ll T) {
	assert(st.mn == nv.mn() && far.mn == nv.mn());
	assert(st.topL() == nv.topL() && st.topR() == nv.topR());
	assert(far.topL() == (nv.topL() == -SlopeTrick::inf ? -SlopeTrick::inf : nv.topL() + T));
	assert(far.topR() == (nv.topR() == SlopeTrick::inf ? SlopeTrick::inf : nv.topR() + T));
	rep(x,-B,B+1) assert(evalST(st, x) == nv.at(x));
}

int main() {
	rep(it,0,60000) {
		SlopeTrick st, far;
		Naive nv;
		int ops = rand() % 15, V2 = rand() % V + 1, shifts = 0;
		ll T = (rand() % 2 ? 1 : -1) * (ll)1e18 - (rand() % 2 ? 0 : 50 + rand() % 50);
		check(st, nv, far, T);
		rep(o,0,ops) {
			int t = rand() % 7;
			ll a = rand() % (2 * V2 + 1) - V2, b = rand() % (2 * V2 + 1) - V2;
			if (t == 0) st.add(a), far.add(a), nv.add(a);
			if (t == 1) st.addR(a), far.addR(a + T), nv.addR(a);
			if (t == 2) st.addL(a), far.addL(a + T), nv.addL(a);
			if (t == 3) st.addAbs(a), far.addAbs(a + T), nv.addAbs(a);
			if (t == 4) st.prefMin(), far.prefMin(), nv.prefMin();
			if (t == 5) st.sufMin(), far.sufMin(), nv.sufMin();
			if (t == 6) {
				if (shifts++ == S) continue;
				a = rand() % (2 * S + 1) - S, b = rand() % (2 * S + 1) - S;
				if (rand() % 3 == 0) b = a; else if (a > b) swap(a, b);
				st.shift(a, b), far.shift(a, b), nv.shift(a, b);
			}
			check(st, nv, far, T);
		}
	}
	// Usage example: min sum |a_i - b_i| with b non-decreasing.
	rep(it,0,20000) {
		int n = rand() % 6, m = rand() % 5 + 1;
		vector<ll> a(n);
		for (ll& x : a) x = rand() % m;
		SlopeTrick f;
		for (ll x : a) f.prefMin(), f.addAbs(x);
		vector<ll> dp(m, 0); // dp[v] = best with b_i = v
		rep(i,0,n) {
			rep(v,1,m) dp[v] = min(dp[v], dp[v - 1]);
			rep(v,0,m) dp[v] += abs(a[i] - v);
		}
		assert(f.mn == *min_element(all(dp)));
	}
	cout << "Tests passed!" << endl;
}
