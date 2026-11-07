#include "../utilities/template.h"
#include "../../content/geometry/PolygonCenter.h"

int main() {
	auto p = polygonCenter(vector<Point<double>>{});
	assert(isnan(p.x) && isnan(p.y));
	auto q = polygonCenter(vector<Point<ll>>{});
	assert(isnan(q.x) && isnan(q.y));
	cout << "Tests passed!" << endl;
}
