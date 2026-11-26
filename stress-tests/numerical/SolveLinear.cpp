#include "../utilities/template.h"
#include "../../content/number-theory/ModInt.h"
#include "../../content/numerical/SolveLinear.h"

template<int P> ll bruteForce(int nmmax) {
	using M = ModInt<P>;
	ll cases = 0;
	rep(n,0,5) rep(m,0,5) {
		if (n * m > nmmax) continue;
		int cells = n * m + n; ll total = 1, xs = 1;
		rep(i,0,cells) total *= P;
		rep(i,0,m) xs *= P;
		for (ll code = 0; code < total; code++) {
			vector<vector<M>> A(n, vector<M>(m));
			vector<M> b(n), x(m, 7);
			ll c = code;
			rep(i,0,n) rep(j,0,m) A[i][j] = c % P, c /= P;
			rep(i,0,n) b[i] = c % P, c /= P;
			auto holds = [&](const vector<M>& v, bool homogeneous) {
				rep(i,0,n) {
					M s = 0;
					rep(j,0,m) s += A[i][j] * v[j];
					if (s != (homogeneous ? M(0) : b[i])) return false;
				}
				return true;
			};
			// Solutions, and which coordinates take a single value.
			ll sols = 0; vector<M> first(m); vector<bool> unique(m, true);
			for (ll xc = 0; xc < xs; xc++) {
				vector<M> v(m); ll t = xc;
				rep(j,0,m) v[j] = t % P, t /= P;
				if (!holds(v, false)) continue;
				if (!sols++) first = v;
				rep(j,0,m) if (v[j] != first[j]) unique[j] = false;
			}
			auto A2 = A; auto b2 = b; vector<vector<M>> ker(3);
			int r = solveLinear(A2, b2, x, &ker);
			A2 = A; b2 = b; vector<M> x1(m);
			assert(solveLinear(A2, b2, x1) == r);
			cases++;
			if (!sols) { assert(r == -1); continue; }
			assert(x1 == x && holds(x, false));
			assert(0 <= r && r <= min(n, m) && sz(ker) == m - r);
			ll want = 1;
			rep(i,0,m-r) want *= P;
			assert(sols == want);
			// ker spans the null space: P^(m-r) distinct combinations, all in it.
			set<vector<int>> seen;
			for (ll kc = 0; kc < want; kc++) {
				vector<M> v(m); ll t = kc;
				for (auto& w : ker) {
					rep(j,0,m) v[j] += M(t % P) * w[j];
					t /= P;
				}
				assert(holds(v, true));
				vector<int> key;
				for (M e : v) key.push_back(e.val());
				seen.insert(key);
			}
			assert((ll)sz(seen) == want);
			rep(j,0,m) {
				bool zero = true;
				for (auto& w : ker) if (w[j] != 0) zero = false;
				assert(zero == unique[j]);
			}
		}
	}
	return cases;
}

typedef vector<double> vd;
double rnd(mt19937& rng) { return uniform_real_distribution<double>(-1, 1)(rng); }

void testDouble() {
	mt19937 rng(12345);
	rep(it,0,2000) {
		int n = (int)(rng() % 8) + 1, m = (int)(rng() % 8) + 1;
		int k = (int)(rng() % (min(n, m) + 1));
		// A = L R has rank k almost surely.
		vector<vd> L(n, vd(k)), R(k, vd(m)), A(n, vd(m));
		for (auto& row : L) for (auto& v : row) v = rnd(rng);
		for (auto& row : R) for (auto& v : row) v = rnd(rng);
		rep(i,0,n) rep(j,0,m) rep(t,0,k) A[i][j] += L[i][t] * R[t][j];
		vd x0(m), b(n), x(m);
		for (auto& v : x0) v = rnd(rng);
		rep(i,0,n) rep(j,0,m) b[i] += A[i][j] * x0[j];
		auto A2 = A; auto b2 = b; vector<vd> ker;
		assert(solveLinear(A2, b2, x, &ker) == k);
		assert(sz(ker) == m - k);
		rep(i,0,n) {
			double s = 0;
			rep(j,0,m) s += A[i][j] * x[j];
			assert(abs(s - b[i]) < 1e-7);
			for (auto& w : ker) {
				s = 0;
				rep(j,0,m) s += A[i][j] * w[j];
				assert(abs(s) < 1e-7);
			}
		}
		if (k == m) rep(j,0,m) assert(abs(x[j] - x0[j]) < 1e-6);
		if (k < n) {
			// b outside the column space: no solution.
			A2 = A; b2 = b;
			vd y(n); // y orthogonal to the columns, from the left null space
			vector<vd> At(m, vd(n)), lker;
			rep(i,0,n) rep(j,0,m) At[j][i] = A[i][j];
			vd zero(m), tmp(n);
			solveLinear(At, zero, tmp, &lker);
			assert(sz(lker) == n - k);
			rep(i,0,n) b2[i] += lker[0][i];
			assert(solveLinear(A2, b2, x) == -1);
		}
	}
	int n = 400;
	vector<vd> A(n, vd(n));
	vd b(n), x(n);
	for (auto& row : A) for (auto& v : row) v = rnd(rng) * 1000;
	for (auto& v : b) v = rnd(rng) * 1000;
	auto A2 = A; auto b2 = b;
	assert(solveLinear(A2, b2, x) == n);
	rep(i,0,n) {
		double s = 0;
		rep(j,0,n) s += A[i][j] * x[j];
		assert(abs(s - b[i]) < 1e-6);
	}
}

// Repeated equations must cancel exactly, and a determined variable
// must have exact zeros in the null space.
void testDoubleExact() {
	{
		vector<vector<double>> A{{7, 21}, {7, 21}, {0, 7}}, ker;
		vector<double> b{10001, 10001, 14}, x(2);
		assert(solveLinear(A, b, x, &ker) == 2 && ker.empty());
		assert(abs(x[1] - 2) < 1e-9 && abs(x[0] - (10001 / 7.0 - 6)) < 1e-9);
	}
	mt19937 rng(5);
	uniform_real_distribution<double> un(-10, 10);
	rep(it,0,20000) {
		// x0 is fixed by row 0; rows 1, 2 are dependent in x1..x3.
		double a = un(rng), c = un(rng), d = un(rng), e = un(rng), f = un(rng);
		vector<vector<double>> A{{a, 0, 0, 0}, {a, c, d, e}, {f * a, f * c, f * d, f * e}}, ker;
		vector<double> b{a, a, f * a}, x(4);
		assert(solveLinear(A, b, x, &ker) == 2 && sz(ker) == 2);
		assert(abs(x[0] - 1) < 1e-9);
		for (auto& w : ker) assert(w[0] == 0);
	}
}

int main() {
	ll cases = bruteForce<2>(12) + bruteForce<3>(8) + bruteForce<5>(4);
	testDouble();
	testDoubleExact();
	cout << "Tests passed! (" << cases << " exhaustive systems)" << endl;
}
