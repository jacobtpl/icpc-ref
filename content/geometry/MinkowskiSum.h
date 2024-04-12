/**
 * Author: Benq
 * Description: Minkowski sum of two convex polygons given in CCW order.
 * example: diameter of hull A is max in A+(-A)
 * Time: O(N)
 * Source: https://cp-algorithms.com/geometry/minkowski.html#implementation
 * Status: tested on 
   * MITPC
   * https://open.kattis.com/problems/roberthood
 */
using vP = vector<P>;
vP minkowski_sum(vP a, vP b) {
	if (sz(a) > sz(b)) swap(a, b);
	if (!sz(a)) return {};
	if (sz(a) == 1) {for (auto &t: b) t += a.front(); return b;}
	rotate(begin(a), min_element(all(a)), end(a));
	rotate(begin(b), min_element(all(b)), end(b));
	a.pb(a[0]), a.pb(a[1]);
	b.pb(b[0]), b.pb(b[1]);
	vP result;
	int i = 0, j = 0;
	while (i < sz(a)-2 || j < sz(b)-2) {
		result.pb(a[i]+b[j]);
		T crs = (a[i+1]-a[i]).cross(b[j+1]-b[j]);
		i += (crs >= 0);
		j += (crs <= 0);
	}
	return result;
}