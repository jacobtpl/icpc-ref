#include "../utilities/template.h"

#include "../../content/numerical/SolveLinearBinary.h"

// Old exhaustive test of a copy of the algorithm using bitset<5>.
namespace old {
const int nmax = 5, mmax = 5, nmmax = 16;

typedef bitset<5> bs;

int solveLinear(vector<bs>& A, vi& b, bs& x, int m) {
	int n = sz(A), rank = 0, br;
	assert(m <= sz(x));
	vi col(m); iota(all(col), 0);
	rep(i,0,n) {
		for (br=i; br<n; ++br) if (A[br].any()) break;
		if (br == n) {
			rep(j,i,n) if(b[j]) return -1;
			break;
		}
		int bc = (int)A[br]._Find_next(i-1);
		swap(A[i], A[br]);
		swap(b[i], b[br]);
		swap(col[i], col[bc]);
		rep(j,0,n) if (A[j][i] != A[j][bc]) {
			A[j].flip(i); A[j].flip(bc);
		}
		rep(j,i+1,n) if (A[j][i]) {
			b[j] ^= b[i];
			A[j] ^= A[i];
		}
		rank++;
	}

	x = bs();
	for (int i = rank; i--;) {
		if (!b[i]) continue;
		x[col[i]] = 1;
		rep(j,0,i) b[j] ^= A[j][i];
	}
	return rank; // (multiple solutions if rank < m)
}

template<class F>
void rec(int i, int j, vector<bs>& A, int m, F f) {
	if (i == sz(A)) {
		f();
	}
	else if (j == m) {
		rec(i+1, 0, A, m, f);
	}
	else {
		rep(v,0,2) {
			A[i][j] = v;
			rec(i, j+1, A, m, f);
		}
	}
}

template<class F>
void rec2(int i, bs& A, int m, F f) {
	if (i == m) f();
	else {
		rep(v,0,2) {
			A[i] = v;
			rec2(i+1, A, m, f);
		}
	}
}

void brute() {
	int ct = 0;
	rep(n,0,nmax+1) rep(m,0,mmax+1) {
		int nm = n*m;
		if (nm > nmmax) continue;
		vector<bs> A(n, bs(m));
		bs b, x, theX;
		vi b2(n);
		rec(0, 0, A, m, [&]() {
			rec2(0, b, n, [&]() {
				int sols = 0;
				rec2(0, x, m, [&]() {
					rep(i,0,n) {
						int v = 0;
						rep(j,0,m) v ^= A[i][j] & x[j];
						if (v != b[i]) return;
					}
					sols++;
					if (sols == 1) theX = x;
				});
				vector<bs> A2 = A;
				bs x2 = x; rep(i,0,n) b2[i] = b[i];
				int r = solveLinear(A2, b2, x2, m);
				if (sols == 0) assert(r == -1);
				else if (sols == 1) assert(r == m);
				else assert(r < m);
				if (sols == 1) assert(x2 == theX);
				ct++;
			});
		});
	}
}
} // namespace old

// ---- tests of the actual header (bitset<1000>) ----
mt19937 rng(99);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

// independent oracle: rank of A and of [A|b] over F_2
int rank2(vector<vi> a) {
	int n = sz(a), m = n ? sz(a[0]) : 0, r = 0;
	rep(c,0,m) {
		int p = r;
		while (p < n && !a[p][c]) p++;
		if (p == n) continue;
		swap(a[p], a[r]);
		rep(i,r+1,n) if (a[i][c]) rep(k,c,m) a[i][k] ^= a[r][k];
		if (++r == n) break;
	}
	return r;
}

void check(const vector<vi>& A, const vi& b, int m) {
	int n = sz(A);
	vector<vi> Ab(n, vi(m + 1));
	vector<bs> A2(n); vi b2 = b; bs x; x.set();
	rep(i,0,n) { rep(j,0,m) { Ab[i][j] = A[i][j]; if (A[i][j]) A2[i][j] = 1; } Ab[i][m] = b[i]; }
	int rk2 = rank2(Ab);
	rep(i,0,n) Ab[i].pop_back();
	int rk = rank2(Ab);
	int r = solveLinear(A2, b2, x, m);
	if (rk != rk2) { assert(r == -1); return; }
	assert(r == rk);
	rep(j,m,sz(x)) assert(!x[j]);
	rep(i,0,n) {
		int s = 0;
		rep(j,0,m) s ^= A[i][j] & (int)x[j];
		assert(s == b[i]);
	}
}

void testRandom(int iters, int maxN, int maxM) {
	rep(it,0,iters) {
		int n = rnd(0, maxN), m = rnd(0, maxM), dens = rnd(1, 4), mode = rnd(0, 2);
		vector<vi> A(n, vi(m)); vi b(n);
		rep(i,0,n) {
			if (mode && i && rnd(0, 1)) { // xor of earlier rows -> rank deficiency
				rep(t,0,2) { int k = rnd(0, i-1); rep(j,0,m) A[i][j] ^= A[k][j]; b[i] ^= b[k]; }
				if (mode == 2) b[i] = rnd(0, 1);
			} else {
				rep(j,0,m) A[i][j] = rnd(0, dens) == 0;
				b[i] = rnd(0, 1);
			}
		}
		check(A, b, m);
	}
}

int main() {
	old::brute();
	{
		vector<bs> A; vi b; bs x; x.set();
		assert(solveLinear(A, b, x, 0) == 0 && x.none());
		assert(solveLinear(A, b, x, 1000) == 0 && x.none());
		A.assign(2, bs()); b = {0, 1};
		assert(solveLinear(A, b, x, 3) == -1);
		A.assign(1, bs()); b = {1}; A[0][999] = 1; // last column only
		assert(solveLinear(A, b, x, 1000) == 1 && x.count() == 1 && x[999]);
	}
	testRandom(200000, 6, 6);
	testRandom(20000, 12, 12);
	testRandom(300, 60, 60);
	rep(it,0,3) { // full width
		int n = it == 0 ? 1000 : rnd(900, 1100), m = 1000;
		vector<vi> A(n, vi(m)); vi b(n);
		rep(i,0,n) { rep(j,0,m) A[i][j] = rnd(0, 1); b[i] = rnd(0, 1); }
		if (it == 2) rep(i,n/2,n) { A[i] = A[i - n/2]; b[i] = b[i - n/2]; rep(j,0,m) A[i][j] ^= A[0][j]; b[i] ^= b[0]; }
		check(A, b, m);
	}
#ifdef SOLVELINEARBINARY_HIGH_BITS
	{ // bits at positions >= m (e.g. b stored as an augmented column) are not ignored:
		// col[bc] is indexed with bc >= m -> out of bounds write
		vector<bs> A(2); vi b = {1, 0}; bs x;
		A[0][2] = 1; A[1][0] = A[1][1] = 1; // x0 + x1 = 0, and row 0 is "0 = 1" plus a stray bit 2
		int r = solveLinear(A, b, x, 2);
		cout << "returned " << r << ", expected -1" << endl;
		assert(r == -1);
	}
#endif
	cout<<"Tests passed!"<<endl;
}
