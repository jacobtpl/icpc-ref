#include "../template.h"
#include "../../content/graph/NetworkSimplex.h"

void print(i128 x) {
	if (x < 0) cout << '-', x = -x;
	string s;
	do s += char('0' + int(x % 10)); while (x /= 10);
	reverse(all(s));
	cout << s << '\n';
}

int main() {
	cin.tie(0)->sync_with_stdio(0);
	int n, m; cin >> n >> m;
	if (!n) { cout << "0\n"; return 0; }
	NetworkSimplex ns; ns.init(n);
	rep(i,0,n) { ll b; cin >> b; ns.B[i] += b; }
	vi id(m);
	vector<ll> lo(m), loop(m);
	i128 extra = 0;
	rep(i,0,m) {
		int a, b; ll u, c; cin >> a >> b >> lo[i] >> u >> c;
		if (a == b) id[i] = -1, loop[i] = c < 0 ? u : lo[i], extra += i128(loop[i]) * c;
		else id[i] = ns.ae(a, b, lo[i], u, c);
	}
	i128 res;
	try { res = ns.solve(); } catch (int) { cout << "infeasible\n"; return 0; }
	print(res + extra);
	rep(i,0,n) cout << ns.P[i] - ns.P[0] << '\n'; // keep |p| small
	rep(i,0,m) cout << (id[i] == -1 ? loop[i] : lo[i] + ns.E[id[i]].flow) << '\n';
}
