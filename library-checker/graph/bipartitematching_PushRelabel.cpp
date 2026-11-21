#include "../template.h"
#include "../../content/graph/PushRelabel.h"

int main() {
	cin.tie(0)->sync_with_stdio(0);
	int L, R, M; cin >> L >> R >> M;
	PushRelabel d(L + R + 2);
	int s = L + R, t = s + 1;
	rep(i,0,M) { int a, b; cin >> a >> b; d.addEdge(a, L + b, 1); }
	rep(i,0,L) d.addEdge(s, i, 1);
	rep(i,0,R) d.addEdge(L + i, t, 1);
	cout << d.calc(s, t) << '\n';
	rep(a,0,L) for (auto& e : d.g[a]) if (e.f > 0 && e.dest < L + R)
		cout << a << ' ' << e.dest - L << '\n';
}
