#include "../utilities/template.h"

#include "../../content/numerical/SimulatedAnnealing.h"

// --- RNG: every (rnd() % n, rnd() % n) pair is reachable, rnd01 in [0,1)
void testRng() {
	for (int n : {1, 2, 64, 256, 500, 1024}) {
		vector<char> seen((size_t)n * n);
		ll draws = 40LL * n * n, hit = 0;
		rep(d,0,draws) {
			int i = int(rnd() % n), j = int(rnd() % n);
			assert(0 <= i && i < n && 0 <= j && j < n);
			if (!seen[i * n + j]) seen[i * n + j] = 1, hit++;
		}
		assert(hit == (ll)n * n);
	}
	int bits[64] = {};
	const int N = 1 << 20;
	double sum = 0;
	rep(d,0,N) {
		uint64_t x = rnd();
		rep(b,0,64) bits[b] += int(x >> b & 1);
		double r = rnd01();
		assert(0 <= r && r < 1);
		sum += r;
	}
	rep(b,0,64) assert(abs(bits[b] - N / 2) < 6 * 512); // 6 sigma
	assert(abs(sum / N - 0.5) < 0.005);
	// any seed works, including 0, and seeding is deterministic
	rs = 0; uint64_t a = rnd(), b = rnd(); assert(a != b);
	rs = 0; assert(rnd() == a && rnd() == b);
	rs = ~0ULL; assert(rnd() != 0 || rnd() != 0);
	rs = 88172645463325252;
}

// --- edge cases of anneal() with counting callbacks
void testEdge() {
	int steps, undos, saves;
	auto reset = [&]() { steps = undos = saves = 0; };
	auto stepC = [&](double v) { return [&, v]() { steps++; return v; }; };
	auto undo = [&]() { undos++; };
	auto save = [&]() { saves++; };

	// sec <= 0: nothing is called, cur is returned
	for (double sec : {0.0, -1.0, -1e300}) {
		reset();
		assert(anneal(5, sec, 10, 1, stepC(0), undo, save) == 5);
		assert(steps == 0 && undos == 0 && saves == 0);
	}
	// constant score: always accepted, never saved
	reset();
	assert(anneal(3, 0.02, 10, 1, stepC(3), undo, save) == 3);
	assert(steps > 0 && undos == 0 && saves == 0);
	// always worse at tiny T: every step undone
	reset();
	assert(anneal(0, 0.02, 1e-9, 1e-9, stepC(1), undo, save) == 0);
	assert(steps > 0 && undos == steps && saves == 0);
	// always worse at huge T: (almost) every step accepted, never saved
	reset();
	assert(anneal(0, 0.02, 1e300, 1e300, stepC(1), undo, save) == 0);
	assert(steps > 0 && undos < steps / 100 && saves == 0);
	// strictly decreasing scores: one save per step, best = last
	reset();
	double x = 1e15;
	double r = anneal(x, 0.02, 1e300, 1e-300,
		[&]() { steps++; return x -= 1; }, undo, save);
	assert(r == x && saves == steps && undos == 0);
	// n=1 neighbourhood of huge scores at T = 1e-300: no NaN, state consistent
	reset();
	int state = 0;
	r = anneal(-1e300, 0.02, 1e-300, 1e-300, [&]() {
		steps++; state ^= 1; return state ? 1e300 : -1e300;
	}, [&]() { undos++; state ^= 1; }, save);
	assert(r == -1e300 && state == 0 && undos == steps && saves == 0);
	// runs for about sec seconds
	auto st = chrono::steady_clock::now();
	anneal(0, 0.05, 1, 1, stepC(0), undo, save);
	double el = chrono::duration<double>(
		chrono::steady_clock::now() - st).count();
	assert(0.05 <= el && el < 0.1);
}

// --- small TSP, 2-opt moves, compare with Held-Karp
void testTsp() {
	int optimal = 0, cases = 0;
	rep(it,0,150) {
		int n = 1 + it % 9;
		vector<vector<ll>> d(n, vector<ll>(n));
		rep(i,0,n) rep(j,i+1,n) d[i][j] = d[j][i] = 1 + rand() % 100;
		// exact
		vector<vector<ll>> dp(1 << n, vector<ll>(n, LLONG_MAX / 2));
		dp[1][0] = 0;
		rep(m,1,1<<n) rep(u,0,n) if (dp[m][u] < LLONG_MAX / 2 && (m >> u & 1))
			rep(v,0,n) if (!(m >> v & 1))
				dp[m | 1 << v][v] = min(dp[m | 1 << v][v], dp[m][u] + d[u][v]);
		ll opt = LLONG_MAX;
		rep(u,0,n) opt = min(opt, dp[(1 << n) - 1][u] + d[u][0]);

		vi p(n), bp;
		iota(all(p), 0);
		random_shuffle(all(p));
		auto len = [&](const vi& q) {
			ll s = 0;
			rep(i,0,n) s += d[q[i]][q[(i + 1) % n]];
			return s;
		};
		int i = 0, j = 0, steps = 0, undos = 0, saves = 0;
		double cur = (double)len(p), lastSaved = cur;
		auto rev = [&](int a, int b) { reverse(p.begin() + a, p.begin() + b + 1); };
		double b = anneal(cur, 0.005, 100, 0.1, [&]() {
			steps++;
			j = int(rnd() % n), i = int(rnd() % (j + 1));
			rev(i, j); return (double)len(p);
		}, [&]() { undos++; rev(i, j); }, [&]() {
			saves++; bp = p;
			assert((double)len(p) < lastSaved);
			lastSaved = (double)len(p);
		});
		assert(steps > 0 && undos <= steps);
		assert(b >= (double)opt);
		if (saves == 0) {
			assert(b == cur);
		} else {
			assert(sz(bp) == n);
			vi s = bp; sort(all(s));
			rep(k,0,n) assert(s[k] == k);
			assert((double)len(bp) == b);
		}
		cases++;
		optimal += (b == (double)opt);
	}
	assert(optimal == cases);
}

// --- maximize by negating: max-cut on a small graph vs brute force
void testMaxCut() {
	rep(it,0,40) {
		int n = 2 + rand() % 11, m = rand() % (n * (n - 1) / 2 + 1);
		vector<array<int, 3>> e;
		rep(k,0,m) {
			int u = rand() % n, v = rand() % n;
			if (u != v) e.push_back({u, v, rand() % 50});
		}
		auto cut = [&](int mask) {
			ll s = 0;
			for (auto& x : e) s += ((mask >> x[0] ^ mask >> x[1]) & 1) * x[2];
			return s;
		};
		ll opt = 0;
		rep(mask,0,1<<n) opt = max(opt, cut(mask));
		int mask = rand() % (1 << n), bm = mask, bit = 0;
		double b = -anneal(-(double)cut(mask), 0.005, 50, 0.1, [&]() {
			bit = int(rnd() % n); mask ^= 1 << bit; return -(double)cut(mask);
		}, [&]() { mask ^= 1 << bit; }, [&]() { bm = mask; });
		assert(b <= (double)opt);
		assert((double)cut(bm) == b);
		assert(b == (double)opt);
	}
}

int main() {
	testRng();
	testEdge();
	testTsp();
	testMaxCut();
	cout << "Tests passed!" << endl;
}
