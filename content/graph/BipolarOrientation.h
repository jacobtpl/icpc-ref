/**
 * Author: Siyong
 * Date: 2024-02-08
 * Source: Tarjan (https://en.wikipedia.org/wiki/Bipolar_orientation#Algorithms)
 * Description: Finds a bipolar orientation of a biconnected graph
 * Time: O(M)
 * Status: https://codeforces.com/contest/1916/submission/245524057
 */

// Warning: Mutates the vector `a`
vector<int> bipolarOrient(vector<vector<int> > a, int s, int t) {
	size_t N = a.size(); // must have s != t, N >= 2
	vector<int> o(N), p(N, -1), d(N, -1), l(N), lk[2];
	lk[0]=lk[1]=vector<int>(N,-1); // lk[0] = prev, lk[1] = next
	a[s].insert(a[s].begin(), t); // can duplicate edge
	int time=0;
	auto f=[&](auto& f, int n) ->void{
		o[time]=n,l[n]=d[n]=time++;
		for(int x:a[n]) if(x!=p[n])
		{
			if(d[x]==-1) {
				p[x]=n, f(f, x); // assert(n==s || l[x]<d[n]);
				ckmin(l[n], l[x]);
			} else ckmin(l[n], d[x]);
		}
	};
	f(f, s);
	auto add=[&](int u, int v, bool b){
		lk[!b][v]=lk[!b][u]; // b true: before, b false: after
		lk[!b][u]=v;
		lk[b][v]=u;
		if(lk[!b][v]!=-1) lk[b][lk[!b][v]]=v;
	};
	add(s, t, 0);
	vector<char> sgn(N, 0);
	sgn[t]=1;
	for(int i=2;i<N;++i) {
		int n=o[i];
		add(p[n], n, sgn[p[n]]=!sgn[o[l[n]]]);
	} // assert(lk[0][s] == -1);
	vector<int> ans;
	for(;s!=-1;s=lk[1][s]) ans.push_back(s);
	return ans;
}
