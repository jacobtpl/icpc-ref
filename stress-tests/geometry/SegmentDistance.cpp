#include "../utilities/template.h"
#include "../../content/geometry/SegmentDistance.h"

int main() {
	const P s(0,0), e(2,0), p(1,1);
	assert(segDist(s,e,p) == 1);
	assert(segDist(P(0,0),P(2,0),P(1,1)) == 1);
	assert(segDist(s,e,P(3,0)) == 1);
	assert(segDist(s,s,P(3,4)) == 5);
	cout << "Tests passed!" << endl;
}
