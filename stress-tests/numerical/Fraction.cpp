#include "../utilities/template.h"
#include "../../content/number-theory/Fraction.h"
#include "../../content/data-structures/Matrix.h"
#include "../../content/numerical/Determinant.h"
#include "../../content/numerical/Gaussian.h"
#include "../../content/numerical/MatrixInverse.h"
#include "../../content/numerical/SolveLinear.h"
#include "../../content/numerical/Tridiagonal.h"
#include "../../content/numerical/PolyInterpolate.h"

using F=Fraction;
using Mat=vector<vector<F>>;

F bruteDet(const Mat& a) {
	int n=sz(a); vi perm(n); iota(all(perm),0);
	F out=0;
	do {
		F term=1;
		rep(i,0,n) {
			term*=a[i][perm[i]];
			rep(j,0,i) if (perm[j]>perm[i]) term=-term;
		}
		out+=term;
	} while (next_permutation(all(perm)));
	return out;
}

vector<F> multiply(const Mat& a, const vector<F>& x) {
	vector<F> b(sz(a));
	rep(i,0,sz(a)) rep(j,0,sz(x)) b[i]+=a[i][j]*x[j];
	return b;
}

int main() {
	mt19937 rng(654123);
	rep(it,0,4000) {
		int n=int(rng()%4)+1, den=int(rng()%5)+1;
		Mat a(n,vector<F>(n));
		for (auto& row : a) for (auto& x : row) x=F(int(rng()%7)-3,den);
		F want=bruteDet(a); Mat copy=a;
		assert(det(copy)==want);
		copy=a; auto [value,rank]=gauss(copy);
		assert(value==want && ((rank==n)==(want!=0)));
		copy=a; assert((matInv(copy)==n)==(want!=0));
		if (want!=0) rep(i,0,n) rep(j,0,n) {
			F sum=0;
			rep(k,0,n) sum+=a[i][k]*copy[k][j];
			assert(sum==(i==j));
		}
		vector<F> x(n), result(n);
		for (F& v : x) v=F(int(rng()%9)-4,3);
		auto b=multiply(a,x), rhs=b; copy=a;
		assert(solveLinear(copy,rhs,result)==rank);
		assert(multiply(a,result)==b);
		if (want!=0) assert(result==x);
	}
	rep(it,0,2000) {
		int n=int(rng()%4)+1, m=int(rng()%5)+1;
		Mat a(n,vector<F>(m)); vector<F> x(m), result(m);
		for (auto& row : a) for (auto& x : row) x=int(rng()%7)-3;
		for (F& v : x) v=F(int(rng()%7)-3,5);
		auto b=multiply(a,x), rhs=b; Mat copy=a;
		int rank=solveLinear(copy,rhs,result);
		assert(rank>=0 && multiply(a,result)==b);
		copy=a; assert(gauss(copy).second==rank);
	}
	Mat tiny={{F(1,1000000000000LL)}};
	assert(gauss(tiny).second==1 && tiny[0][0]==1);
	tiny={{F(1,1000000000000LL)}};
	assert(matInv(tiny)==1 && tiny[0][0]==1000000000000LL);
	tiny={{F(1,1000000000000LL)}};
	vector<F> rhs={1}, result(1);
	assert(solveLinear(tiny,rhs,result)==1 && result[0]==1000000000000LL);
	Mat inconsistent={{1,1},{2,2}};
	rhs={1,3}; result.resize(2);
	assert(solveLinear(inconsistent,rhs,result)==-1);
	rep(it,0,1000) {
		int n=int(rng()%5)+1;
		vector<F> diag(n,4), upper(n-1,F(1,2)), lower(n-1,F(-1,3));
		vector<F> want(n), b(n);
		for (F& v : want) v=F(int(rng()%9)-4,5);
		rep(i,0,n) {
			b[i]=diag[i]*want[i];
			if (i) b[i]+=lower[i-1]*want[i-1];
			if (i+1<n) b[i]+=upper[i]*want[i+1];
		}
		assert(tridiagonal(diag,upper,lower,b)==want);
		vector<F> x(n), y(n), coefficients(n);
		rep(i,0,n) coefficients[i]=F(int(rng()%7)-3,3);
		rep(i,0,n) {
			x[i]=i;
			for (int j=n; j--;) y[i]=y[i]*x[i]+coefficients[j];
		}
		assert(interpolate(x,y,n)==coefficients);
	}
	assert(tridiagonal(vector<F>{0,0},vector<F>{1},vector<F>{1},
		vector<F>{2,3})==vector<F>({3,2}));
	Matrix<F,2> a;
	a.d={{{F(1,2),F(1,2)},{F(1,2),F(1,2)}}};
	assert((a^100).d==a.d);
	assert((a^0)*vector<F>({2,3})==vector<F>({2,3}));
	// Floating-point instantiations retain their tolerance policies.
	vector<vector<double>> d={{1e-15}};
	assert(gauss(d).second==0);
	d={{1e-15}}; assert(matInv(d)==0);
	d={{1e-15}}; vector<double> b={1}, x(1);
	assert(solveLinear(d,b,x)==-1);
	cout << "Fraction algebra tests passed (7000 random cases)\n";
}
