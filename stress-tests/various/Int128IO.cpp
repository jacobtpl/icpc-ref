#include "../utilities/template.h"
#include "../../content/various/Int128IO.h"

using I = __int128_t;
using U = __uint128_t;

void check(bool ok) { if (!ok) abort(); }

int main(int argc, char**) {
	if (argc > 1) {
		string s;
		while (getline(cin,s)) {
			I x = 42;
			if (parse128(s,x)) cout << str128(x) << '\n';
			else { check(x == 42); cout << "invalid\n"; }
		}
		return 0;
	}
	const U sign = U(1) << 127;
	const I hi = I(sign-1), lo = -hi-1;
	check(str128(lo) == "-170141183460469231731687303715884105728");
	check(str128(hi) == "170141183460469231731687303715884105727");
	I x = 0;
	check(parse128(str128(lo),x) && x == lo);
	check(parse128(str128(hi),x) && x == hi);
	for (string s : {"", "+", "-", " 1", "1 ", "1.0", "--1", "1a",
		"170141183460469231731687303715884105728",
		"-170141183460469231731687303715884105729"}) {
		x = 42;
		check(!parse128(s,x) && x == 42);
	}
	check(parse128("+00042",x) && x == 42);
	check(parse128("-00000",x) && x == 0 && str128(x) == "0");
	for (ll n = -50000; n <= 50000; ++n) {
		check(str128(I(n)) == to_string(n));
		check(parse128(to_string(n),x) && x == n);
	}
	mt19937_64 rng(128);
	for (int i = 0; i < 100000; ++i) {
		U bits = (U(rng()) << 64) | rng();
		I n = I(bits & (sign-1));
		if (bits & sign) n = -n-1;
		check(parse128(str128(n),x) && x == n);
	}
	cout << "Int128 I/O tests passed (100000 wide, 100001 small)\n";
}
