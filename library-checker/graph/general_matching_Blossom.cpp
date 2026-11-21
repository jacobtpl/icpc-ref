#include "../template.h"
#include "../../content/graph/GeneralMatchingBlossom.h"

int main() {
	cin.tie(0)->sync_with_stdio(0);
	int n, m; cin >> n >> m;
	MaxMatching mm; mm.init(n);
	rep(i,0,m) { int a, b; cin >> a >> b; mm.ae(a + 1, b + 1); }
	cout << mm.solve() << '\n';
	rep(i,1,n+1) if (mm.mate[i] > i) cout << i - 1 << ' ' << mm.mate[i] - 1 << '\n';
}
