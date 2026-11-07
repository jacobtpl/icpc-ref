#include "../utilities/template.h"
#include "../utilities/random.h"

#include "../../content/various/SMAWK.h"

typedef vector<vector<ll>> Mat;

vi brute(const Mat& A) {
	vi ans;
	for (auto& row : A) ans.push_back(int(min_element(all(row)) - row.begin()));
	return ans;
}
vi bruteMax(const Mat& A) { // rightmost maximum
	vi ans;
	for (auto& row : A) {
		int b = 0;
		rep(j,0,sz(row)) if (row[j] >= row[b]) b = j;
		ans.push_back(b);
	}
	return ans;
}
bool totallyMonotone(const Mat& A, int n, int m) {
	rep(i,0,n) rep(i2,i+1,n) rep(j,0,m) rep(j2,j+1,m)
		if (A[i][j] > A[i][j2] && !(A[i2][j] > A[i2][j2])) return false;
	return true;
}
bool monge(const Mat& A, int n, int m) {
	rep(i,0,n-1) rep(j,0,m-1)
		if (A[i][j] + A[i+1][j+1] > A[i][j+1] + A[i+1][j]) return false;
	return true;
}
ll calls;
vi run(const Mat& A, int n, int m) {
	return smawk(n, m, [&](int i, int j) {
		assert(0 <= i && i < n && 0 <= j && j < m);
		calls++;
		return A[i][j];
	});
}
void check(const Mat& A, int n, int m) {
	assert(totallyMonotone(A, n, m));
	vi got = run(A, n, m), want = brute(A);
	if (got != want) {
		cout << "n=" << n << " m=" << m << endl;
		for (auto& row : A) { for (ll x : row) cout << x << ' '; cout << endl; }
		cout << "want:"; for (int x : want) cout << ' ' << x; cout << endl;
		cout << "got: "; for (int x : got) cout << ' ' << x; cout << endl;
		abort();
	}
}
Mat randMonge(int n, int m, int maxD, int maxOff) {
	Mat A(n, vector<ll>(m, 0));
	rep(i,0,n) rep(j,0,m) {
		ll d = randRange(maxD + 1);
		A[i][j] = -d + (i ? A[i-1][j] : 0) + (j ? A[i][j-1] : 0) - (i && j ? A[i-1][j-1] : 0);
	}
	vector<ll> ro(n), co(m);
	for (ll& x : ro) x = randRange(-maxOff, maxOff + 1);
	for (ll& x : co) x = randRange(-maxOff, maxOff + 1);
	rep(i,0,n) rep(j,0,m) A[i][j] += ro[i] + co[j];
	assert(monge(A, n, m));
	return A;
}

int main() {
	// Edge cases
	assert(smawk(0, 0, [](int, int) { return 0; }).empty());
	assert(smawk(0, 5, [](int, int) { return 0; }).empty());
	assert(smawk(1, 1, [](int, int) { return 0; }) == vi{0});
	assert(smawk(3, 1, [](int, int) { return 0; }) == vi({0, 0, 0}));
	assert(smawk(4, 4, [](int, int) { return 0; }) == vi({0, 0, 0, 0}));
	assert(smawk(1, 6, [](int, int j) { return abs(j - 3); }) == vi{3});
	assert(smawk(1, 6, [](int, int j) { return j == 4 ? 0 : 1; }) == vi{4});
	assert(smawk(6, 6, [](int i, int j) { return (i-j)*(i-j); }) == vi({0,1,2,3,4,5}));
	{ // extreme values, no arithmetic in the routine
		const ll B = ll(4e18);
		vi got = smawk(3, 3, [&](int i, int j) { return i == j ? -B : B; });
		assert(got == vi({0, 1, 2}));
		vector<ll> a = {-ll(2e9), 0, ll(2e9)}, b = {-ll(2e9), 0, ll(2e9)};
		got = smawk(3, 3, [&](int i, int j) { return (a[i]-b[j])*(a[i]-b[j]); });
		assert(got == vi({0, 1, 2}));
	}
	{ // doubles
		vector<double> a = {0.1, 0.5, 0.9};
		vi got = smawk(3, 5, [&](int i, int j) { return (a[i] - j / 4.0) * (a[i] - j / 4.0); });
		assert(got == vi({0, 2, 4}));
	}

	// Random Monge matrices, values small so ties are common
	rep(it,0,60000) {
		int n = randRange(1, 10), m = randRange(1, 10);
		Mat A = randMonge(n, m, randRange(3), randRange(4));
		check(A, n, m);
		// row maxima recipes from the header
		Mat B(n, vector<ll>(m));
		rep(i,0,n) rep(j,0,m) B[i][j] = -A[i][m-1-j];
		vi got = run(B, n, m), want = bruteMax(A);
		rep(i,0,n) got[i] = m - 1 - got[i];
		assert(got == want);
		// A2 = -A is inverse Monge, so negating it is valid and gives A2's leftmost maxima
		Mat A2(n, vector<ll>(m));
		rep(i,0,n) rep(j,0,m) A2[i][j] = -A[i][j];
		vi gotMax = smawk(n, m, [&](int i, int j) { return -A2[i][j]; });
		rep(i,0,n) {
			int b = 0;
			rep(j,0,m) if (A2[i][j] > A2[i][b]) b = j;
			assert(gotMax[i] == b);
		}
	}

	// Random totally monotone but non-Monge matrices (rejection sampling)
	int tm = 0;
	rep(it,0,300000) {
		int n = randRange(1, 6), m = randRange(1, 6);
		Mat A(n, vector<ll>(m));
		rep(i,0,n) rep(j,0,m) A[i][j] = randRange(4);
		if (!totallyMonotone(A, n, m)) continue;
		tm++;
		check(A, n, m);
	}
	assert(tm > 10000);

	// Large: (a_i - b_j)^2 + v_j with sorted a, b; compare to D&C monotone minima, check call count
	{
		int n = 300000, m = 300000;
		vector<ll> a(n), b(m), v(m);
		for (ll& x : a) x = randRange(1, 1000000000);
		for (ll& x : b) x = randRange(1, 1000000000);
		for (ll& x : v) x = randRange(0, 1000000);
		sort(all(a)); sort(all(b));
		auto f = [&](int i, int j) { return (a[i]-b[j])*(a[i]-b[j]) + v[j]; };
		calls = 0;
		vi got = smawk(n, m, [&](int i, int j) { calls++; return f(i, j); });
		assert(calls <= 8 * (n + m));
		vi want(n);
		auto rec = [&](auto& self, int L, int R, int lo, int hi) -> void {
			if (L >= R) return;
			int mid = (L + R) / 2, best = lo;
			rep(j,lo+1,hi) if (f(mid, j) < f(mid, best)) best = j;
			want[mid] = best;
			self(self, L, mid, lo, best + 1);
			self(self, mid + 1, R, best, hi);
		};
		rec(rec, 0, n, 0, m);
		assert(got == want);
		// n >> m and m >> n
		got = smawk(n, 100, f);
		rep(i,0,n) { int best = 0; rep(j,1,100) if (f(i, j) < f(i, best)) best = j; assert(got[i] == best); }
		got = smawk(50, m, f);
		rep(i,0,50) { int best = 0; rep(j,1,m) if (f(i, j) < f(i, best)) best = j; assert(got[i] == best); }
	}
	cout << "Tests passed!" << endl;
}
