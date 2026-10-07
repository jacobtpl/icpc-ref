#include "../utilities/template.h"

#define main main2
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#include "../../content/strings/Hashing-codeforces.h"
#pragma GCC diagnostic pop
#undef main

// Independent oracle: the same polynomial mod 1e9+7, 1e9+9 and 2^32,
// packed the way operator ull() packs it.
ull oracle(const string& s, int a, int b) {
	const ll M1 = 1000000007, M2 = 1000000009;
	ll h1 = 0, h2 = 0; unsigned h3 = 0;
	rep(i,a,b) {
		ll c = s[i];
		h1 = ((h1 * C + c) % M1 + M1) % M1;
		h2 = ((h2 * C + c) % M2 + M2) % M2;
		h3 = h3 * (unsigned)C + (unsigned)c;
	}
	return (ull)h1 ^ ((ull)h2 ^ (ull)h3 << 21) << 21;
}

void testOracle(int lo, int hi) { // characters drawn from [lo, hi]
	mt19937 rng(555);
	rep(it,0,3000) {
		int n = (int)(rng() % 50);
		string s(n, 'a');
		bool few = rng() % 2;
		for (char& c : s) c = (char)(few ? lo + rng() % 2 : lo + rng() % (hi - lo + 1));
		HashInterval h(s);
		H whole = hashString(s);
		assert((ull)whole == oracle(s, 0, n));
		vector<pair<H, string>> subs;
		rep(i,0,n+1) rep(j,i,n+1) {
			H x = h.hashInterval(i, j);
			assert((ull)x == oracle(s, i, j));
			if (sz(subs) < 60) subs.emplace_back(x, s.substr(i, j - i));
		}
		// operator== / operator< as they would be used in a solution (map keys, sorting)
		for (auto& x : subs) for (auto& y : subs) {
			assert((x.first == y.first) == (x.second == y.second));
			assert((x.first < y.first) == ((ull)x.first < (ull)y.first));
		}
		map<H, string> m;
		for (auto& x : subs) assert(m.emplace(x.first, x.second).first->second == x.second);
		rep(le,1,n+2) {
			auto ve = getHashes(s, le);
			assert(sz(ve) == max(0, n - le + 1));
			rep(i,0,sz(ve)) assert(ve[i] == h.hashInterval(i, i + le));
		}
	}
}

void testThueMorse() {
	// Thue-Morse strings collide mod 2^32 (the reason for the two extra primes)
	rep(k,0,15) {
		string a(1 << k, 'a'), b = a;
		rep(i,0,1<<k) a[i] = (char)('a' + __builtin_parity(i)), b[i] = (char)('b' - __builtin_parity(i));
		assert(!(hashString(a) == hashString(b)));
	}
}

#include <sys/time.h>
int main() {
	// fixed bases first (deterministic), then the time-based one from the header
	for (int c : {1000003, 999999, 31337, 257}) {
		C = c;
		testOracle(1, 127);
		testThueMorse();
#ifdef TEST_HIGH_CHARS
		// Undocumented precondition: characters must be non-negative
		// (with signed char, bytes >= 128 give non-canonical residues).
		testOracle(-128, 127);
#endif
	}
	timeval tp;
	gettimeofday(&tp, 0);
	C = (int)tp.tv_usec;
	assert((ull)(H(1)*2+1-3) == 0);

	rep(it,0,10000) {
		int n = rand() % 10;
		int alpha = rand() % 10 + 1;
		string s;
		rep(i,0,n) s += (char)('a' + rand() % alpha);
		HashInterval hi(s);
		set<string> strs;
		set<ull> hashes;

		// HashInterval
		rep(i,0,n+1) rep(j,i,n+1) {
			string sub = s.substr(i, j - i);
			ull hash = (ull) hashString(sub);
			assert((ull) hi.hashInterval(i, j) == hash);
			hashes.insert(hash);
			strs.insert(sub);
		}

		// getHashes
		rep(le,1,n+1) {
			auto ve = getHashes(s, le);
			assert(sz(ve) == n-le+1);
			rep(i,0,n-le+1) {
				assert((ull) ve[i] == (ull) hi.hashInterval(i, i + le));
			}
		}

		// No collisions
		assert(sz(strs) == sz(hashes));
	}
	cout<<"Tests passed!"<<endl;
}
