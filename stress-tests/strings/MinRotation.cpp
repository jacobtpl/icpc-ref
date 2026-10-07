#include "../utilities/template.h"

#include "../../content/strings/MinRotation.h"

int min_rotation2(string& v) {
	int n = sz(v);
	string w = v; w.insert(w.end(), all(v));
	int j = 0;
	rep(i,1,n) {
		if (vi(w.begin() + i, w.begin() + i + n) <
			vi(w.begin() + j, w.begin() + j + n)) j = i;
	}
	return j;
}

void testPerf() {
	string s;
	rep(i,0,10000000)
		s += (char)(rand()%400000 < 2);
	cout << minRotation(s) << endl;
}

int main() {
	rep(it,0,1000000) {
		int n = rand() % 10;
		string v;
		rep(i,0,n) v += (char)(rand() % 3);
		int r = minRotation(v);
		int r2 = min_rotation2(v);
		assert(r == r2);
		rotate(v.begin(), v.begin() + r, v.end());
		assert(minRotation(v) == 0);
		assert(min_rotation2(v) == 0);
	}
	assert(minRotation("") == 0);
	// longer strings: periodic, all-equal, near-periodic, high bytes
	mt19937 rng(31337);
	rep(it,0,100000) {
		int n = 1 + (int)(rng() % 60), kind = (int)(rng() % 4);
		int alpha = kind == 0 ? 256 : 1 + (int)(rng() % 3);
		int p = 1 + (int)(rng() % 6);
		string v(n, 0);
		rep(i,0,n) v[i] = (char)('a' + rng() % (unsigned)alpha);
		if (kind >= 2) rep(i,p,n) v[i] = v[i-p];
		if (kind == 3) v[rng() % (unsigned)n]++;
		int r = minRotation(v), best = 0;
		string w = v + v;
		// first index of the smallest rotation, chars compared as char
		rep(i,1,n) if (lexicographical_compare(w.begin()+i, w.begin()+i+n,
			w.begin()+best, w.begin()+best+n)) best = i;
		assert(r == best);
	}
	cout<<"Tests passed!"<<endl;
}
