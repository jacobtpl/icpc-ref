using ll = long long;

int const MOD = 998244353;

struct mint
{
public:
	int v;

	explicit operator int() {return v;}

	mint(): v(0) {}
	template<int> mint(int _v): v(_v)
	{
		if(v < -MOD || MOD <= v)
			v %= MOD;
		if(v < 0) v += MOD;
	}
	template<typename T> mint(T _v)
	{
		v = _v % MOD;
		if(v < 0) v += MOD;
	}

	mint& operator+= (mint const& o) {if((v+=o.v)>=MOD) v-=MOD; return *this;}
	mint& operator-= (mint const& o) {if((v-=o.v)<0) v+=MOD; return *this;}
	mint& operator*= (mint const& o) {v=(ll)v*o.v%MOD; return *this;}

	friend mint operator+ (mint a, mint const& b) {return a+=b;}
	friend mint operator- (mint a, mint const& b) {return a-=b;}
	friend mint operator* (mint a, mint const& b) {return a*=b;}

	template<typename T> friend mint pow(mint a, T b)
	{
		mint r(1);
		for(;b;b>>=1, a*=a)
			if(b&1)
				r *= a;
		return r;
	}
};
