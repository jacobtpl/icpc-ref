#include "../template.h"
#include "../../content/graph/Dinic.h"

int main() {
	cin.tie(0)->sync_with_stdio(0);
	int L, R, M; cin >> L >> R >> M;
	Dinic d(L + R + 2);
	int s = L + R, t = s + 1;
	vector<pii> ed(M);
	for (auto& [a, b] : ed) cin >> a >> b, d.addEdge(a, L + b, 1);
	rep(i,0,L) d.addEdge(s, i, 1);
	rep(i,0,R) d.addEdge(L + i, t, 1);
	cout << d.calc(s, t) << '\n';
	rep(i,0,M) if (d.flow(i)) cout << ed[i].first << ' ' << ed[i].second << '\n';
}
