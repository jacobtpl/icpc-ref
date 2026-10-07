#include "../utilities/template.h"
#define pb push_back // from content/contest/template.cpp

#include "../../content/combinatorial/DeBruijnSeq.h"

ll ipow(ll k, int n) { ll r = 1; while (n--) r *= k; return r; }

// every length-n string over [0,k) appears exactly once cyclically
void check(int k, int n) {
	vi s = deBruijnSeq(k, n);
	ll len = ipow(k, n);
	assert(sz(s) == len);
	for (int x : s) assert(0 <= x && x < k);
	vector<bool> seen(len);
	ll h = 0, top = len / k;
	rep(i,0,n) h = h * k + s[i % len];
	rep(i,0,len) {
		assert(!seen[h]); seen[h] = 1;
		h = (h - s[i] * top) * k + s[(i + n) % len];
	}
}

// lexicographically smallest de Bruijn sequence by exhaustive search
int K, N, LEN; vi cur, bestSeq; vector<bool> used;
bool dfs(int i) {
	if (i == LEN) {
		// check wrap-around windows
		vector<bool> u(LEN);
		rep(s,0,LEN) {
			int h = 0;
			rep(j,0,N) h = h * K + cur[(s + j) % LEN];
			if (u[h]) return 0;
			u[h] = 1;
		}
		bestSeq = cur; return 1;
	}
	rep(c,0,K) {
		cur[i] = c;
		int h = -1;
		if (i >= N - 1) {
			h = 0;
			rep(j,0,N) h = h * K + cur[i - N + 1 + j];
			if (used[h]) continue;
			used[h] = 1;
		}
		if (dfs(i + 1)) return 1;
		if (h >= 0) used[h] = 0;
	}
	return 0;
}
vi brute(int k, int n) {
	K = k; N = n; LEN = (int)ipow(k, n);
	cur.assign(LEN, 0); used.assign(LEN, 0);
	assert(dfs(0));
	return bestSeq;
}

int main() {
	// all (k, n) with k^n <= 2e5
	rep(k,1,41) rep(n,1,19) {
		if (k > 1 && ipow(k, n) > 200000) break;
		check(k, n);
	}
	check(1, 30); check(200, 1); check(300, 2); check(2, 20);
	// FKM yields the lexicographically smallest sequence
	vector<pii> small = {{1,1},{1,3},{2,1},{2,2},{2,3},{2,4},{2,5},
		{3,1},{3,2},{3,3},{4,1},{4,2},{5,2},{6,1}};
	for (auto [k, n] : small) assert(deBruijnSeq(k, n) == brute(k, n));
	cout << "Tests passed!" << endl;
}
