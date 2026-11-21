#include "../template.h"
#include "../../content/graph/hopcroftKarp.h"

int main() {
	cin.tie(0)->sync_with_stdio(0);
	int L, R, M; cin >> L >> R >> M;
	vector<vi> g(L);
	rep(i,0,M) { int a, b; cin >> a >> b; g[a].push_back(b); }
	vi btoa(R, -1);
	cout << hopcroftKarp(g, btoa) << '\n';
	rep(b,0,R) if (btoa[b] != -1) cout << btoa[b] << ' ' << b << '\n';
}
