#include "../utilities/template.h"

#include "../../content/various/CDQ.h"

typedef array<int, 3> P;

vi brute(const vector<P>& p) {
	int n = sz(p);
	vi res(n);
	rep(i,0,n) rep(j,0,n) if (i != j &&
		p[j][0] <= p[i][0] && p[j][1] <= p[i][1] && p[j][2] <= p[i][2])
			res[i]++;
	return res;
}

void check(const vector<P>& p) {
	vi got = cdq(p), want = brute(p);
	if (got != want) {
		cout << "n = " << sz(p) << endl;
		for (auto& a : p) cout << a[0] << ' ' << a[1] << ' ' << a[2] << endl;
		rep(i,0,sz(p)) cout << got[i] << ' ';
		cout << endl;
		rep(i,0,sz(p)) cout << want[i] << ' ';
		cout << endl;
		abort();
	}
}

int rnd(int lo, int hi) { return lo + rand() % (hi - lo + 1); }

int main() {
	check({});
	check({{5, 5, 5}});
	check({{1, 2, 3}, {1, 2, 3}});
	check({{INT_MIN, INT_MIN, INT_MIN}, {INT_MAX, INT_MAX, INT_MAX},
	       {INT_MIN, INT_MAX, INT_MIN}, {INT_MAX, INT_MIN, INT_MAX},
	       {INT_MIN, INT_MIN, INT_MIN}, {0, 0, 0}});
	const int ranges[][2] = {{0, 1}, {0, 3}, {-5, 5}, {-1000000000, 1000000000}};
	rep(it,0,200000) {
		int n = rnd(0, 12), shape = rnd(0, 5);
		if (it % 50 == 0) n = rnd(100, 200);
		auto r = ranges[rnd(0, 3)];
		vector<P> p(n);
		rep(i,0,n) {
			P a;
			rep(k,0,3) a[k] = rnd(r[0], r[1]);
			if (shape == 1) a = {i, i, i};
			if (shape == 2) a = {i, n - i, i};
			if (shape == 3) a[0] = a[1] = a[2];
			if (shape == 4 && i && rand() % 2) a = p[rnd(0, i - 1)];
			if (shape == 5) a[0] &= 1, a[2] &= 1;
			p[i] = a;
		}
		check(p);
	}
	vector<P> big(200000);
	for (auto& a : big) a = {0, 0, 0};
	vi res = cdq(big);
	for (int x : res) assert(x == sz(big) - 1);
	rep(i,0,sz(big)) big[i] = {i, i, i};
	res = cdq(big);
	rep(i,0,sz(big)) assert(res[i] == i);
	cout << "Tests passed!" << endl;
}
