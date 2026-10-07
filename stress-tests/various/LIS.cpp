#include "../utilities/template.h"

#include "../../content/various/LIS.h"

template<class I> vi lisWeak(const vector<I>& S) {
	if (S.empty()) return {};
	vi prev(sz(S));
	typedef pair<I, int> p;
	vector<p> res;
	rep(i,0,sz(S)) {
		// 0 -> i for longest non-decreasing subsequence
		auto it = lower_bound(all(res), p{S[i], i});
		if (it == res.end()) res.emplace_back(), it = res.end()-1;
		*it = {S[i], i};
		prev[i] = it == res.begin() ? 0 : (it-1)->second;
	}
	int L = sz(res), cur = res.back().second;
	vi ans(L);
	while (L--) ans[L] = cur, cur = prev[cur];
	return ans;
}

// O(N^2) DP for the length of the longest (weakly) increasing subsequence
template<class I> int lisLen(const vector<I>& v, bool weak) {
	int n = sz(v), best = 0;
	vi dp(n, 1);
	rep(i,0,n) {
		rep(j,0,i) if (weak ? !(v[i] < v[j]) : v[j] < v[i]) dp[i] = max(dp[i], dp[j] + 1);
		best = max(best, dp[i]);
	}
	return best;
}

template<class I> void check(const vector<I>& v) {
	rep(weak,0,2) {
		vi inds = weak ? lisWeak(v) : lis(v);
		assert(sz(inds) == lisLen(v, weak));
		rep(i,0,sz(inds)) assert(0 <= inds[i] && inds[i] < sz(v));
		rep(i,0,sz(inds)-1) {
			assert(inds[i] < inds[i+1]);
			if (weak) assert(!(v[inds[i+1]] < v[inds[i]]));
			else assert(v[inds[i]] < v[inds[i+1]]);
		}
	}
}

void testLarge() {
	mt19937_64 rng(5);
	auto rnd = [&](ll lo, ll hi) { return uniform_int_distribution<ll>(lo, hi)(rng); };
	check(vi{});
	check(vi{5});
	check(vi{INT_MIN, INT_MAX, INT_MIN, INT_MAX});
	check(vector<ll>{LLONG_MAX, LLONG_MIN, LLONG_MAX, 0, LLONG_MIN});
	rep(it,0,3000) {
		int n = (int)rnd(0, 300), mode = it % 6;
		ll lo = mode == 0 ? -3 : mode == 1 ? LLONG_MIN : mode == 2 ? LLONG_MAX - 5 : mode == 3 ? LLONG_MIN : -n;
		ll hi = mode == 0 ? 3 : mode == 1 ? LLONG_MAX : mode == 2 ? LLONG_MAX : mode == 3 ? LLONG_MIN + 5 : n;
		vector<ll> v(n);
		for (auto& x : v) x = rnd(lo, hi);
		if (mode == 5) sort(all(v));
		if (mode == 4 && it % 12 == 4) sort(v.rbegin(), v.rend());
		check(v);
		vi vi_(n); vector<double> vd(n); vector<pii> vp(n); vector<string> vs(n);
		rep(i,0,n) {
			vi_[i] = (int)(v[i] % 1000), vd[i] = (double)(v[i] % 7) / 2;
			vp[i] = {(int)(v[i] % 3), (int)(v[i] % 5)};
			vs[i] = string(1 + (size_t)(v[i] & 1), char('a' + (v[i] & 6)));
		}
		check(vi_); check(vd); check(vp); check(vs);
	}
}

int main() {
	testLarge();
	rep(weak,0,2) {
		auto lt = [weak](int a, int b) { return weak ? a <= b : a < b; };
		rep(it,0,1000000) {
			int n = rand() % 7;
			vi v(n);
			for(auto &x: v) x = rand() % 4;
			vi inds = weak ? lisWeak(v) : lis(v);
			rep(i,0,sz(inds)-1) {
				assert(lt(v[inds[i]], v[inds[i+1]]));
			}
			rep(bi,0,(1 << n)) {
				int si = (int)bitset<32>(bi).count();
				if (si <= sz(inds)) continue;
				int prev = INT_MIN;
				rep(i,0,n) if (bi & (1 << i)) {
					if (!lt(prev, v[i])) goto next;
					prev = v[i];
				}
				cout << "exists lis of size " << si << " but found only " << sz(inds) << endl;
				for(auto &x: v) cout << x << ' ';
				cout << endl;
				abort();
	next:;
			}
		}
	}
	cout << "Tests passed!" << endl;
}
