/**
 * Author: Siyong
 * Date: April 6, 2024
 * Description: Modular class
 * Status: Tested many many times
 * Time: Faster than kactl mod. Slower than using ll directly
 */
int const MOD = 998244353;

ll euclid(ll a, ll b, ll &x, ll &y) {
	if (!b) return x = 1, y = 0, a;
	ll d = euclid(b, a % b, y, x);
	return y -= a/b * x, d;
}

struct mint {
	int v;
	explicit operator int() {return v;}
	mint(): v(0) {}
	mint(auto z) {
		z %= MOD;
		if (z < 0) z += MOD;
		v = z;
	}
	friend mint invert(mint a) {
		ll x, y, g = euclid(a.v, MOD, x, y);
		assert(g == 1); return mint(x);
	}
	mint& operator+= (mint const& o) {if((v+=o.v)>=MOD) v-=MOD; return *this;}
	mint& operator-= (mint const& o) {if((v-=o.v)<0) v+=MOD; return *this;}
	mint& operator*= (mint const& o) {v=(ll)v*o.v%MOD; return *this;}
	mint& operator/= (mint const& o) {return *this *= invert(o);}

	friend mint operator+ (mint a, mint const& b) {return a+=b;}
	friend mint operator- (mint a, mint const& b) {return a-=b;}
	friend mint operator* (mint a, mint const& b) {return a*=b;}
	friend mint operator/ (mint const& a, mint const& b) {return a*invert(b);}

	friend mint pow(mint a, auto b) {
		mint r(1);
		for(;b;b>>=1, a*=a)
			if(b&1)
				r *= a;
		return r;
	}
};
