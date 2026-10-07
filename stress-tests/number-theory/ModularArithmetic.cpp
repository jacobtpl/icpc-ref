#include "../utilities/template.h"

#include "../../content/number-theory/ModularArithmetic.h"

int main() {
	rep(a,0,mod) rep(b,1,mod) {
		Mod ma(a);
		Mod mb(b);
		Mod mc = ma / mb;
		assert((mc * mb).x == a);
	}
	Mod a = 2;
	ll cur=1;
	rep(i, 0, 17) {
		assert((a ^ i).x == cur);
		cur = (cur * 2) % mod;
		// cout << i << ": " << (a ^ i).x << endl;
	}
	// all operations, exhaustively, for reduced operands
	rep(x,0,mod) rep(y,0,mod) {
		Mod mx(x), my(y);
		assert((mx + my).x == (x + y) % mod);
		assert((mx - my).x == ((x - y) % mod + mod) % mod);
		assert((mx * my).x == x * y % mod);
		if (y) {
			ll q = (mx / my).x;
			assert(0 <= q && q < mod && q * y % mod == x);
			ll iv = mx.invert(my).x;
			assert(0 <= iv && iv < mod && iv * y % mod == 1);
		}
	}
	rep(x,0,mod) {
		ll c = 1;
		rep(e,0,200) {
			assert((Mod(x) ^ e).x == c);
			c = c * x % mod;
		}
		for (ll e : {(ll)1e9, (ll)1e18, LLONG_MAX, (ll)mod - 1, (ll)mod * mod}) {
			ll r = 1, b = x;
			for (ll f = e; f; f /= 2, b = b * b % mod) if (f & 1) r = r * b % mod;
			assert((Mod(x) ^ e).x == r);
		}
	}
	// operands that are not already in [0, mod): negative and large values
	auto norm = [](ll v) { return (v % mod + mod) % mod; };
	vector<ll> vals;
	rep(v,-3*(int)mod,3*(int)mod+1) vals.push_back(v);
	for (ll v : {(ll)1e9, (ll)-1e9, (ll)1e18, (ll)-1e18, (ll)INT_MAX, (ll)INT_MIN, (ll)3e18, (ll)-3e18})
		vals.push_back(v);
	for (ll x : vals) for (ll y : vals) {
		Mod mx(x), my(y);
		assert((mx + my).x == norm(norm(x) + norm(y)));
		assert((mx - my).x == norm(norm(x) - norm(y)));
		assert((mx * my).x == norm(norm(x) * norm(y)));
		if (norm(y)) assert(((mx / my) * my).x == norm(x));
	}
	for (ll x : vals) rep(e,0,20) {
		ll c = 1;
		rep(i,0,e) c = c * norm(x) % mod;
		assert((Mod(x) ^ e).x == c);
	}
	cout<<"Tests passed!"<<endl;
}
