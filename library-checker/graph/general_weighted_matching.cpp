#include "../template.h"
#include "../../content/graph/GeneralWeightedMatching.h"

WeightedMatch<501> wm;

int main() {
	cin.tie(0)->sync_with_stdio(0);
	int n, m; cin >> n >> m;
	wm.init(n);
	rep(i,0,m) { int a, b, w; cin >> a >> b >> w; wm.ae(a + 1, b + 1, w); }
	ll w = wm.calc().first;
	int cnt = 0;
	rep(i,1,n+1) cnt += wm.match[i] > i;
	cout << cnt << ' ' << w << '\n';
	rep(i,1,n+1) if (wm.match[i] > i) cout << i - 1 << ' ' << wm.match[i] - 1 << '\n';
}
