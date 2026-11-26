#include "../utilities/template.h"

#include "../../content/numerical/SolveLinearBinary.h"

const int nmax = 5, mmax = 5, nmmax = 16;

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

int main() {
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
	cout<<"Tests passed!"<<endl;
}
