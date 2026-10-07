#include "../utilities/template.h"

#include "../../content/strings/Hashing.h"

// Independent oracle: polynomial hash mod 2^64-1, reduced with %.
const ull MOD = ~0ULL;
ull oracle(const string& s, int a, int b) {
	__uint128_t h = 0, c = C.x % MOD;
	rep(i,a,b) {
		// same value as the implicit char -> ull conversion in H(str[i])
		ull v = (ull)s[i] % MOD;
		h = (h * c + v) % MOD;
	}
	return (ull)h;
}

void testOracle() {
	mt19937_64 rng(987);
	// arithmetic near the wrap-around points
	vector<ull> vals = {0, 1, 2, MOD, MOD - 1, MOD - 2, 1ULL << 63, (1ULL << 63) - 1, (1ULL << 32) - 1, 1ULL << 32};
	rep(i,0,200) vals.push_back(rng());
	for (ull a : vals) for (ull b : vals) {
		__uint128_t x = a % MOD, y = b % MOD;
		assert((H(a) + H(b)).get() % MOD == (ull)((x + y) % MOD));
		assert((H(a) - H(b)).get() % MOD == (ull)((x + MOD - y) % MOD));
		assert((H(a) * H(b)).get() % MOD == (ull)(x * y % MOD));
		assert((H(a) == H(b)) == (x == y));
		assert((H(a) < H(b)) == (H(a).get() < H(b).get()));
	}
	rep(it,0,3000) {
		int n = (int)(rng() % 60);
		string s(n, 'a');
		int mode = (int)(rng() % 3);
		for (char& c : s) c = (char)(mode == 0 ? 'a' + rng() % 2 : mode == 1 ? rng() % 256 : 127 + rng() % 3);
		HashInterval hi(s);
		assert(hashString(s).get() % MOD == oracle(s, 0, n));
		rep(i,0,n+1) rep(j,i,n+1)
			assert(hi.hashInterval(i, j).get() % MOD == oracle(s, i, j));
		rep(le,1,n+2) {
			auto ve = getHashes(s, le);
			assert(sz(ve) == max(0, n - le + 1));
			rep(i,0,sz(ve)) assert(ve[i] == hi.hashInterval(i, i + le));
		}
	}
	// Thue-Morse: the two complementary strings collide mod 2^64 from length 2^10 on
	rep(k,0,15) {
		string a(1 << k, 'a'), b = a;
		rep(i,0,1<<k) a[i] = (char)('a' + __builtin_parity(i)), b[i] = (char)('b' - __builtin_parity(i));
		assert(!(hashString(a) == hashString(b)));
		ull x = 0, y = 0;
		for (char c : a) x = x * C.x + (ull)c;
		for (char c : b) y = y * C.x + (ull)c;
		assert((x == y) == (k >= 10)); // sanity check of the claim in the header
	}
	// no collisions among all substrings of a long random binary string of a fixed length
	{
		string s(200000, 'a');
		for (char& c : s) c = (char)('a' + rng() % 2);
		for (int le : {40, 1000}) {
			auto ve = getHashes(s, le);
			map<ull, int> seen;
			rep(i,0,sz(ve)) {
				auto it = seen.emplace(ve[i].get(), i).first;
				assert(s.compare(it->second, le, s, i, le) == 0);
			}
		}
	}
}

#include <sys/time.h>
int main() {
	assert((H(1)*2+1-3).get() == 0);
	testOracle();

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
			ull hash = hashString(sub).get();
			assert(hi.hashInterval(i, j).get() == hash);
			hashes.insert(hash);
			strs.insert(sub);
		}

		// getHashes
		rep(le,1,n+1) {
			auto ve = getHashes(s, le);
			assert(sz(ve) == n-le+1);
			rep(i,0,n-le+1) {
				assert(ve[i].get() == hi.hashInterval(i, i + le).get());
			}
		}

		// No collisions
		assert(sz(strs) == sz(hashes));
	}
	cout<<"Tests passed!"<<endl;
}
