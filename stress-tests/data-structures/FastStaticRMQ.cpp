#include "../utilities/template.h"

#include "../../content/data-structures/FastStaticRMQ.h"

// query(l, r) is the minimum over the INCLUSIVE range [l, r].
template<class T, class F>
void exhaustive(int N, F gen) {
	vector<T> v(N);
	rep(i,0,N) v[i] = gen();
	RMQ<T> rmq(v);
	rep(i,0,N) {
		T m = v[i];
		rep(j,i,N) {
			m = min(m, v[j]);
			assert(rmq.query(i, j) == m);
		}
	}
}

template<class T, class F>
void randomQueries(int N, int Q, F gen, mt19937& rng) {
	vector<T> v(N);
	rep(i,0,N) v[i] = gen();
	RMQ<T> rmq(v);
	vector<vector<T>> sp(1, v);
	for (int k = 1; (1 << k) <= N; k++) {
		sp.emplace_back(N - (1 << k) + 1);
		rep(i,0,sz(sp[k]))
			sp[k][i] = min(sp[k-1][i], sp[k-1][i + (1 << (k-1))]);
	}
	rep(it,0,Q) {
		int l = (int)(rng() % N), r = (int)(rng() % N);
		if (it % 4 == 0) r = min(N - 1, l + (int)(rng() % 70));
		if (it % 64 == 1) l = 0;
		if (it % 64 == 2) r = N - 1;
		if (l > r) swap(l, r);
		int k = 31 - __builtin_clz(r - l + 1);
		assert(rmq.query(l, r) == min(sp[k][l], sp[k][r - (1 << k) + 1]));
	}
}

int main() {
	mt19937 rng(5);
	{ vi e; RMQ<int> rmq(e); } // empty input must not crash
	rep(N,1,200) {
		exhaustive<int>(N, [&]() { return (int)(rng() % 3); }); // many ties
		exhaustive<int>(N, [&]() { return (int)(rng() % 1000) - 500; });
		exhaustive<int>(N, [&]() { return (int)rng(); }); // full int range
		exhaustive<ll>(N, [&]() { return rng() % 2 ? LLONG_MAX - (ll)(rng() % 5) : LLONG_MIN + (ll)(rng() % 5); });
		exhaustive<double>(N, [&]() { return (double)(rng() % 100) / 7 - 5; });
		int c = N; exhaustive<int>(N, [&]() { return c--; }); // decreasing
		c = 0; exhaustive<int>(N, [&]() { return c++; }); // increasing
		exhaustive<int>(N, [&]() { return 7; });
	}
	rep(it,0,300) {
		int N = (int)(rng() % 3000) + 1;
		exhaustive<int>(N > 600 ? N % 600 + 1 : N, [&]() { return (int)(rng() % 50); });
		randomQueries<int>(N, 3000, [&]() { return (int)(rng() % 20); }, rng);
		randomQueries<ll>(N, 3000, [&]() { return (ll)(((uint64_t)rng() << 32) ^ rng()); }, rng);
	}
	for (int N : {29, 30, 31, 59, 60, 61, 899, 900, 901, 1 << 15, 30 << 10, (30 << 10) - 1, 200000, 1000000}) {
		randomQueries<int>(N, 200000, [&]() { return (int)rng(); }, rng);
		randomQueries<int>(N, 200000, [&]() { return (int)(rng() % 2); }, rng);
		int c = N; randomQueries<int>(N, 200000, [&]() { return c--; }, rng);
		c = 0; randomQueries<int>(N, 200000, [&]() { return c++; }, rng);
	}
	cout<<"Tests passed!"<<endl;
}
