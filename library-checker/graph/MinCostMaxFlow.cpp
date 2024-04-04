#include <bits/stdc++.h>
#include <bits/extc++.h> /// include-line, keep-include

using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

bool ckmax(auto &a, auto const& b) {return b>a?a=b,1:0;}
bool ckmin(auto &a, auto const& b) {return b<a?a=b,1:0;}

struct MCMF {
	typedef ll C; typedef ll F; typedef __int128_t R;
	C const INFC = numeric_limits<C>::max() / 4;
	F const INFF = numeric_limits<F>::max() / 4;

	struct Edge {int n; F flow; C cost; size_t rev;};
	int N;
	vector<vector<Edge>> ed;
	vector<C> dist, pi;
	vector<F> amt;
	vector<size_t> par;

	MCMF(int N) :
		N(N), ed(N), dist(N), pi(N), amt(N), par(N) {}

	void addEdge(int u, int v, F f, C c) {
		ed[u].emplace_back(v, -f, c, ed[v].size());
		ed[v].emplace_back(u, 0, -c, ed[u].size()-1);
	}

	void path(int s) {
		fill(all(amt), 0); amt[s] = INFF;
		fill(all(dist), INFC); dist[s] = 0;

		__gnu_pbds::priority_queue<pair<C, int>> q;
		vector<decltype(q)::point_iterator> its(N);
		q.push({0, s});

		while(!q.empty()) {
			s = q.top().second; q.pop();
			C di = dist[s] + pi[s];
			for(auto [n, f, c, rev]: ed[s])
				if(f < 0 && ckmin(dist[n], di + c - pi[n])) {
					par[n] = rev; amt[n] = min(amt[s], -f);
					if(its[n] == q.end()) its[n] = q.push({-dist[n], n});
					else q.modify(its[n], {-dist[n], n});
				}
		}
		rep(i,0,N) pi[i] = amt[i] ? pi[i] + dist[i] : INFC;
	}

	pair<F, R> maxflow(int s, int t) {
		F totflow = 0; R totcost = 0;
		while (path(s), amt[t]) {
			F fl = amt[t]; totflow += fl;
			for(int n = t; n != s;) {
				auto &[p, f, c, rev] = ed[n][par[n]];
				f -= fl; ed[p][rev].flow += fl; n = p;
				totcost -= (R) fl * c;
			}
		}
		return {totflow, totcost};
	}

	void setpi() {
		fill(all(pi), 0);
		int it = N, ch = 1;
		while(ch-- && it--)
			rep(i,0,N) if(pi[i] != INFC)
				for(auto [to, f, c, _]: ed[i])
					if(f < 0 && ckmin(pi[to], pi[i] + c)) ch = 1;
		assert(it >= 0);
	}
};

void fail() {
	printf("infeasible\n");
	exit(0);
}

void print(__int128_t v) {
	if(v < 0) v = -v, cout << '-';
	vector<int> d;
	for(;v > 0;v /= 10) d.push_back(v % 10);
	if(d.empty()) d.resize(1);
	reverse(all(d));
	for(int x: d) cout << x;
	cout << '\n';
}
int main() {
	cin.tie(0)->sync_with_stdio(0);
	cin.exceptions(cin.failbit);

	int N, M; cin >> N >> M;
	MCMF mf(N + 2);

	vector<ll> B(N);
	for(int i = 0;i < N; ++i)
		cin >> B[i];
	vector<int> S(M), T(M), L(M), U(M), C(M), id(M);

	__int128_t min_cost = 0;
	for(int i =0;i < M; ++i)
	{
		cin >> S[i] >> T[i] >> L[i] >> U[i] >> C[i];
		if(C[i] < 0) {
			min_cost += (ll)C[i] * U[i];
			B[S[i]] -= U[i];
			B[T[i]] += U[i];
			id[i] = mf.ed[T[i]].size();
			mf.addEdge(T[i], S[i], U[i] - L[i], -C[i]);
		}
		else {
			B[S[i]] -= L[i];
			B[T[i]] += L[i];
			min_cost += (ll) C[i] * L[i];
			id[i] = mf.ed[S[i]].size();
			mf.addEdge(S[i], T[i], U[i] - L[i], C[i]);
		}
	}

	ll net_pos_flow = 0;
	ll net_flow = 0;
	for(int i = 0;i < N; ++i)
	{
		if(B[i] > 0) {
			net_pos_flow += B[i];
			mf.addEdge(N, i, B[i], 0);
		} else
			mf.addEdge(i, N+1, -B[i], 0);
		net_flow += B[i];
	}
	if(net_flow != 0) fail();

	auto [max_flow, extra_cost] = mf.maxflow(N, N+1);
	if(max_flow != net_pos_flow) fail();

	print(min_cost + extra_cost);
	mf.setpi();
	for(int i = 0;i < N; ++i) cout << mf.pi[i] << '\n';
	for(int i = 0;i < M; ++i) {
		if(C[i] < 0) {
			ll flow = mf.ed[T[i]][id[i]].flow;
			ll orig_flow = U[i] - L[i];
			cout << U[i] - (flow + orig_flow) << '\n';
		} else {
			ll flow = mf.ed[S[i]][id[i]].flow;
			ll orig_flow = U[i] - L[i];
			cout << L[i] + (flow + orig_flow) << '\n';
		}
	}
	return 0;
}

