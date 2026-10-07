#include "../utilities/template.h"

// KnuthDP.h is documentation only; this implements the DP exactly as
// described there and checks it against the O(N^3) recurrence, for
// random cost functions f satisfying the stated sufficient criteria.

mt19937 rng(777);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

const ll INF = LLONG_MAX / 4;

// a[i][i] = a[i][i+1] = 0 (no valid k); p is the minimal optimal k.
void knuth(int n, const vector<vector<ll>>& f, vector<vector<ll>>& a) {
	vector<vi> p(n+1, vi(n+1));
	a.assign(n+1, vector<ll>(n+1, 0));
	rep(i,0,n) p[i][i+1] = i;
	rep(len,2,n+1) rep(i,0,n+1-len) {
		int j = i + len;
		a[i][j] = INF;
		rep(k,max(i+1,p[i][j-1]),min(j-1,p[i+1][j])+1) {
			ll v = a[i][k] + a[k][j];
			if (v < a[i][j]) a[i][j] = v, p[i][j] = k;
		}
		assert(a[i][j] < INF);
		a[i][j] += f[i][j];
	}
}

void brute(int n, const vector<vector<ll>>& f, vector<vector<ll>>& a, vector<vi>& p) {
	a.assign(n+1, vector<ll>(n+1, 0));
	p.assign(n+1, vi(n+1));
	rep(len,2,n+1) rep(i,0,n+1-len) {
		int j = i + len;
		a[i][j] = INF;
		rep(k,i+1,j) {
			ll v = a[i][k] + a[k][j];
			if (v < a[i][j]) a[i][j] = v, p[i][j] = k;
		}
		a[i][j] += f[i][j];
	}
}

// Random f with f(b,c) <= f(a,d) and
// f(a,c) + f(b,d) <= f(a,d) + f(b,c) for a <= b <= c <= d.
vector<vector<ll>> genF(int n, int type) {
	vector<vector<ll>> f(n+1, vector<ll>(n+1));
	if (type == 0) { // sum of weights (optimal BST / file merging)
		vector<ll> pre(n+1);
		rep(i,0,n) pre[i+1] = pre[i] + rnd(0, type == 0 ? 20 : 0);
		rep(i,0,n+1) rep(j,i,n+1) f[i][j] = pre[j] - pre[i];
	} else if (type == 1) { // squared sum
		vector<ll> pre(n+1);
		rep(i,0,n) pre[i+1] = pre[i] + rnd(0, 5);
		rep(i,0,n+1) rep(j,i,n+1) f[i][j] = (pre[j] - pre[i]) * (pre[j] - pre[i]);
	} else { // sum of nonnegative weights over cells of sub-rectangles
		vector<vector<ll>> w(n+2, vector<ll>(n+2));
		// f(i,j) = sum of w[x][y] for i <= x <= y <= j
		rep(x,0,n+1) rep(y,x,n+1) w[x][y] = rnd(0, 3) ? 0 : rnd(0, 9);
		rep(i,0,n+1) rep(j,i,n+1) {
			ll s = 0;
			rep(x,i,j+1) rep(y,x,j+1) s += w[x][y];
			f[i][j] = s;
		}
	}
	rep(a,0,n+1) rep(b,a,n+1) rep(c,b,n+1) rep(d,c,n+1) {
		assert(f[b][c] <= f[a][d]);
		assert(f[a][c] + f[b][d] <= f[a][d] + f[b][c]);
	}
	return f;
}

int main() {
	rep(it,0,60000) {
		int n = it < 50000 ? rnd(1, 9) : rnd(1, 30);
		auto f = genF(n, it % 3);
		vector<vector<ll>> a, b; vector<vi> p;
		knuth(n, f, a);
		brute(n, f, b, p);
		assert(a == b);
		// the claimed monotonicity of the minimal optimal k
		rep(i,0,n+1) rep(j,i+3,n+1) {
			assert(p[i][j-1] <= p[i][j]);
			assert(p[i][j] <= p[i+1][j]);
		}
	}
	// O(N^2) claim: count inner iterations at N = 2000
	{
		int n = 2000;
		vector<ll> pre(n+1);
		rep(i,0,n) pre[i+1] = pre[i] + rnd(0, 1000);
		vector<vi> p(n+1, vi(n+1));
		vector<vector<ll>> a(n+1, vector<ll>(n+1, 0));
		rep(i,0,n) p[i][i+1] = i;
		ll ops = 0;
		rep(len,2,n+1) rep(i,0,n+1-len) {
			int j = i + len;
			a[i][j] = INF;
			rep(k,max(i+1,p[i][j-1]),min(j-1,p[i+1][j])+1) {
				ll v = a[i][k] + a[k][j]; ops++;
				if (v < a[i][j]) a[i][j] = v, p[i][j] = k;
			}
			a[i][j] += pre[j] - pre[i];
		}
		assert(ops <= 4LL * n * n);
	}
	cout << "Tests passed!" << endl;
}
