#include "../template.h"
#include "../../content/graph/WeightedMatching.h"

int main() {
	cin.tie(0)->sync_with_stdio(0);
	int n; cin >> n;
	vector<vector<ll>> a(n, vector<ll>(n));
	for (auto& r : a) for (auto& x : r) cin >> x;
	auto [cost, match] = hungarian(a);
	cout << cost << '\n';
	rep(i,0,n) cout << match[i] << " \n"[i == n - 1];
}
