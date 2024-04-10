// Builds Suffix Tree, then Builds Suffix Arr
#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

bool ckmax(auto &a, auto const& b) {return b>a?a=b,1:0;}
bool ckmin(auto &a, auto const& b) {return b<a?a=b,1:0;}

struct SA {
	vector<map<char, int> > adj;
	vi link, dis;
	int N;
	SA(): adj(1), link(1, -1), dis(1, 0), N(1) {}
	int new_node(int v=-1) {
		if(v == -1)
			adj.emplace_back(), link.emplace_back(), dis.emplace_back();
		else
			adj.push_back(adj[v]), link.push_back(link[v]), dis.push_back(dis[v]);
		return N++;
	}
	int go(int p, int c) {
		auto it = adj[p].find(c);
		if(dis[it->second] == dis[p] + 1)
			return it->second;
		else {
			int q = it->second, n = new_node(q);
			dis[n] = dis[p] + 1, link[q] = n;
			for(;p != -1 && (it = adj[p].find(c))->second == q;p = link[p])
				it->second = n;
			return n;
		}
	}
	int append(int p, char c) {
		auto it = adj[p].find(c);
		if(it != adj[p].end())
			return go(p, c);
		int n = new_node();
		dis[n] = dis[p] + 1;
		for(;p != -1 && adj[p].find(c) == adj[p].end();p = link[p])
			adj[p].insert({c, n});
		if(p == -1)
			link[n] = 0;
		else
			link[n] = go(p, c);
		return n;
	}
	int add(string const &s) {
		int n = 0;
		for(char c: s)
			n = append(n, c);
		return n;
	}
};

int main() {
	cin.tie(0)->sync_with_stdio(0);
	cin.exceptions(cin.failbit);

	SA sa;
	string s;
	cin >> s;
	s += '$';
	int N = s.size();

	vector<int> sloc(N+1);
	for(int i = 0;i < N; ++i)
		sloc[i] = sa.append(sloc[i-1], s[N-i-1]);

	vector<map<char, tuple<int, int, int> > > treeAdj(sa.N); // node, [s[l]...s[r])
	vector<char> vis(sa.N);
	vis[0] = 1;

	for(int i = 1;i <= N; ++i) { // tree edges store last occurrence
	// for(int i = N;i >= 1; --i) { // tree edges store first occurrence
		int n = sloc[i], p;

		int r = N;
		for(;;n=p) {
			if(vis[n]) break;
			vis[n] = 1;
			p = sa.link[n];

			int len = sa.dis[n] - sa.dis[p];
			treeAdj[p].insert({s[r-len], {n, r-len, r}});
			r -= len;
		}
	}
	
	vector<int> remap(sa.N, -1);
	for(int i = 1;i <= N; ++i)
		remap[sloc[i]] = N - i;
	vi ans;

	auto dfs = [&](auto &&self, int n, int d = 0) ->void{
		if(treeAdj[n].empty())
			ans.push_back(N - d);
		for(auto [c, info]: treeAdj[n]){
			auto [x, l, r] = info;
			// printf("%d -> %c[%d %d %d]: %s\n", n, c, x, l, r, s.substr(l, r-l).c_str());
			self(self, x, d + r-l);
		}
	};
	dfs(dfs, 0);

	for(int i = 0;i < ans.size(); ++i)
		printf("%d%c", ans[i], " \n"[i+1==ans.size()]);
	return 0;
}

