/**
 * Author: Siyong
 * Date: April 6, 2024
 * Description: Boilerplate centroid decomp code
 * Status: Tested by Jacob on Xenia and Tree
 * Time: O(N\log N)
 */

struct Centroid {
	vector<vi> const& adj;
	int N;
	vector<int> s, rem;
	vector<vi> links;
	vector<pair<int, int> > par; // <parent centroid, index>
	int root;
	int dfs(int n, int p=-1) {
		s[n]=1;
		for(int x: adj[n]) if(x!=p && !rem[x])
			s[n]+=dfs(x,n);
		return s[n];
	}
	int find(int n, int t, int p=-1) {
		for(int k=1;k--;)
			for(int x:adj[n]) if(x!=p && !rem[x] && s[x]*2>t)
			{p=n,n=x,k=1; break;}
		return n;
	}
	int cent(int start=0) {
		int c = find(start, dfs(start));
		// Do stuff with c. Just remember to check both (x != p && !rem[x])
		rem[c]=1;
		for(int x:adj[c]) if(!rem[x]) {
			int v = cent(x);
			par[v] = {c, sz(links[c])};
			links[c].push_back(v);
		}
		return c;
	}
	Centroid(vector<vi> const& adj): adj(adj), N(adj.size()), s(N), rem(N), links(N), par(N,{-1,-1}), root(cent()) {}
};
