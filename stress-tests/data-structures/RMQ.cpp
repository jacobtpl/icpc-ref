#include "../utilities/template.h"

#include "../../content/data-structures/RMQ.h"

int main(int argc, char**) {
	srand(2);
	rep(N,0,100) {
		vi v(N);
		rep(i,0,N) v[i] = i;
		random_shuffle(all(v));
		RMQ<int> rmq(v);
		rep(i,0,N) rep(j,i+1,N+1) {
			int m = rmq.query(i,j);
			int n = 1 << 29;
			rep(k,i,j) n = min(n, v[k]);
			assert(n == m);
		}
	}
	// duplicates, negatives, extreme values, other value types
	mt19937_64 rng(5);
	rep(it,0,3000) {
		int N = (int)(rng() % 70) + 1, mode = (int)(rng() % 4);
		vector<ll> v(N);
		for (auto& x : v)
			x = mode == 0 ? (ll)(rng() % 3) - 1 : mode == 1 ? (ll)rng() :
				mode == 2 ? (rng() % 2 ? LLONG_MIN : LLONG_MAX) : -(ll)(rng() % 1000);
		RMQ<ll> rmq(v);
		rep(i,0,N) {
			ll m = LLONG_MAX;
			rep(j,i+1,N+1) {
				m = min(m, v[j-1]);
				assert(rmq.query(i,j) == m);
			}
		}
	}
	{
		RMQ<int> one(vi{-7}), none(vi{});
		assert(one.query(0, 1) == -7);
		assert(sz(none.jmp) == 1 && none.jmp[0].empty());
		vector<pair<double, string>> ps = {{1.5, "b"}, {1.5, "a"}, {-2, "z"}, {3, ""}};
		RMQ<pair<double, string>> pr(ps);
		assert(pr.query(0, 2).second == "a" && pr.query(0, 4).second == "z" && pr.query(3, 4).first == 3);
	}
	for (int N : {127, 128, 129, 1000, 4096, 100000}) { // powers of two +-1, large
		vi v(N);
		for (int& x : v) x = (int)rng();
		RMQ<int> rmq(v);
		rep(it,0,2000) {
			int a = (int)(rng() % N), b = (int)(rng() % N);
			if (a > b) swap(a, b);
			if (it < 4) a = it & 1 ? 0 : a, b = it & 2 ? N - 1 : b;
			assert(rmq.query(a, b + 1) == *min_element(v.begin() + a, v.begin() + b + 1));
		}
	}
	if (argc > 1) { // benchmark mode: ./a.out bench
		for (int N : {1000000, 5000000}) {
			vi v(N);
			for (int& x : v) x = (int)rng();
			auto t0 = chrono::steady_clock::now();
			RMQ<int> rmq(v);
			double tb = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
			vector<pii> qs(N);
			for (auto& q : qs) {
				q.first = (int)(rng() % N), q.second = (int)(rng() % N);
				if (q.first > q.second) swap(q.first, q.second);
				q.second++;
			}
			t0 = chrono::steady_clock::now();
			unsigned h = 0;
			for (auto& q : qs) h = h * 31 + rmq.query(q.first, q.second);
			cerr << "N=Q=" << N << ": build " << tb << " s, queries " << chrono::duration<double>(
				chrono::steady_clock::now() - t0).count() << " s (" << h << ")\n";
		}
	}
	cout<<"Tests passed!"<<endl;
}
