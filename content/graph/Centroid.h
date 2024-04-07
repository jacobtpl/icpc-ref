struct Centroid {
	vector<vi> const& adj;
	int N;
	vector<int> s, rem;
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
	void cent(int start=0) {
		int c = find(start, dfs(start));
		// Do stuff with c. Just remember to check both (x != p && !rem[x])
		rem[c]=1;
		for(int x:adj[c]) if(!rem[x])
			cent(x);
	}
	Centroid(vector<vi> const& adj): adj(adj), N(adj.size()), s(N), rem(N) {
		cent();
	}
};
