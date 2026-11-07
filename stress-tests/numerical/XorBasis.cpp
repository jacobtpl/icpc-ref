#include "../utilities/template.h"

#include "../../content/numerical/XorBasis.h"

typedef unsigned long long ull;

mt19937_64 rng(1337);
ull randVal(int mode) {
	ull x = rng();
	switch (mode) {
	case 0: return x % 8;
	case 1: return x % 256;
	case 2: return x;
	case 3: return 1ULL << (x % 64);
	case 4: return ~0ULL - x % 4;
	case 5: return (ull)(-(ll)(x % 100));
	case 6: return (x >> 60 << 60) | (x & 15);
	default: return x & (x >> 7) & (x >> 13);
	}
}

vector<ull> allXors(const vector<ull>& v) {
	set<ull> s = {0};
	for (ull x : v) {
		set<ull> t = s;
		for (ull y : s) t.insert(y ^ x);
		s = t;
	}
	return vector<ull>(all(s));
}

void check(XorBasis& B, const vector<ull>& v) {
	vector<ull> rep = allXors(v);
	assert((1ULL << B.r) == rep.size());
	int cnt = 0;
	rep(i,0,64) if (B.b[i]) {
		cnt++;
		assert(63 - __builtin_clzll(B.b[i]) == i);
		assert(binary_search(all(rep), B.b[i]));
	}
	assert(cnt == B.r);
	rep(k,0,sz(rep)) assert(B.kth((ull)k) == rep[k]);
	assert(B.maxXor() == rep.back());
	assert(B.minXor(0) == 0 && B.has(0));
	rep(it,0,8) {
		ull x = it < 4 ? rep[rng() % rep.size()] : randVal((int)(rng() % 8));
		assert(B.has(x) == binary_search(all(rep), x));
		ull mn = ~0ULL, mx = 0;
		for (ull y : rep) mn = min(mn, x ^ y), mx = max(mx, x ^ y);
		assert(B.minXor(x) == mn);
		assert(B.maxXor(x) == mx);
	}
}

int main() {
	{
		XorBasis B;
		assert(B.r == 0 && B.has(0) && !B.has(1));
		assert(B.kth(0) == 0 && B.maxXor() == 0 && B.minXor(7) == 7 && B.maxXor(7) == 7);
		assert(!B.add(0) && B.r == 0);
		assert(B.add(5) && B.add(3) && !B.add(6) && B.r == 2);
		assert(B.has(6) && B.maxXor() == 6 && B.kth(2) == 5);
		XorBasis E; E.merge(B); assert(E.r == 2 && E.maxXor() == 6);
		B.merge(B); assert(B.r == 2 && B.kth(3) == 6);
	}
	{
		XorBasis B;
		rep(i,0,64) assert(B.add(1ULL << i));
		assert(B.r == 64 && !B.add(rng()) && B.maxXor() == ~0ULL);
		rep(it,0,1000) {
			ull k = rng();
			assert(B.kth(k) == k && B.has(k) && B.minXor(k) == 0 && B.maxXor(k) == ~0ULL);
		}
		XorBasis C;
		rep(i,0,64) C.add(rng() | 1ULL << 63);
		rep(i,0,200) C.add(rng());
		assert(C.r == 64 && C.kth(12345) == 12345);
	}
	rep(mode,0,8) rep(it,0,1500) {
		int n = (int)(rng() % 9);
		vector<ull> v;
		XorBasis B;
		rep(i,0,n) {
			ull x = randVal(mode);
			if (rng() % 5 == 0 && !v.empty()) x = v[rng() % v.size()];
			if (rng() % 7 == 0) x = 0;
			vector<ull> before = allXors(v);
			bool added = B.add(x);
			assert(added == !binary_search(all(before), x));
			v.push_back(x);
			check(B, v);
		}
		int m = (int)(rng() % 5);
		vector<ull> w;
		XorBasis C;
		rep(i,0,m) { w.push_back(randVal(mode)); C.add(w.back()); }
		XorBasis Ccopy = C;
		B.merge(C);
		assert(memcmp(C.b, Ccopy.b, sizeof C.b) == 0 && C.r == Ccopy.r);
		for (ull x : w) v.push_back(x);
		check(B, v);
	}
	cout<<"Tests passed!"<<endl;
}
