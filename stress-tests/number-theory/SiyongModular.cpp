#include "../utilities/template.h"
#include <sys/wait.h>
#include <unistd.h>

#include "../../content/number-theory/SiyongModular.h"

template<class T> bool ckmax(T& a, T b) {
	return b > a ? (a = b, true) : false;
}
namespace gaussian {
#include "../../content/numerical/Gaussian.h"
}

using i128 = __int128_t;
using u128 = __uint128_t;

static_assert(sizeof(mint) == sizeof(int));
static_assert(mint::mod() == MOD);
static_assert(mint::raw(5).val() == 5);
static_assert(int(-mint(1)) == MOD-1);
static_assert(+mint(-1) == MOD-1);
static_assert(0 == mint(0) && mint(1) != 0);
static_assert(pow(mint(2), 10) == 1024);
static_assert(!is_constructible_v<ModInt<7>, ModInt<11>>);
static_assert(!is_constructible_v<mint, double>);

template<int M> int reduce(i128 x) {
	x %= M;
	if (x < 0) x += M;
	return int(x);
}
template<int M> int power(int a, u128 b) {
	int r = 1 % M;
	for (; b; b /= 2, a = reduce<M>(i128(a)*a))
		if (b & 1) r = reduce<M>(i128(r)*a);
	return r;
}
template<int M> int inverse(int a) {
	if constexpr (M <= 17) {
		for (int b = 0; b < M; ++b)
			if (a*b % M == 1 % M) return b;
		abort();
	} else return power<M>(a, M-2); // large test moduli are prime
}
template<int M> void checkPair(ll x, ll y) {
	using Z = ModInt<M>;
	const Z a(x), b(y);
	int av = reduce<M>(x), bv = reduce<M>(y);
	assert(a.v == av && b.v == bv);
	assert(int(a) == av && a.val() == av);
	assert((a == b) == (av == bv));
	assert((a != b) == (av != bv));
	assert((a == 0) == (av == 0) && (0 == a) == (av == 0));
	assert((a != 0) == (av != 0) && (0 != a) == (av != 0));
	assert((a+b).v == reduce<M>(i128(av)+bv));
	assert((a-b).v == reduce<M>(i128(av)-bv));
	assert((a*b).v == reduce<M>(i128(av)*bv));
	assert((-a).v == reduce<M>(-i128(av)) && +a == a);
	assert(x+b == a+y && x+b == a+b);
	assert(x-b == a-y && x-b == a-b);
	assert(x*b == a*y && x*b == a*b);
	Z c = a;
	assert(&(c += b) == &c && c == a+b);
	c = a; assert(&(c -= b) == &c && c == a-b);
	c = a; assert(&(c *= b) == &c && c == a*b);
	c = a; c += c; assert(c == a+a);
	c = a; c -= c; assert(c == 0);
	c = a; c *= c; assert(c == a*a);
	c = a; assert(c++ == a && c == a+1);
	assert(c-- == a+1 && c == a);
	assert(&++c == &c && c == a+1);
	assert(&--c == &c && c == a);
	--c; assert(c == a-1);
	if (gcd(bv, M) == 1) {
		int ib = inverse<M>(bv);
		assert(invert(b).v == ib);
		assert((a/b).v == reduce<M>(i128(av)*ib));
		assert(x/b == a/y && x/b == a/b);
		c = a; assert(&(c /= b) == &c && c == a/b);
		c = b; c /= c; assert(c == 1);
	}
	assert(a.v == av && b.v == bv);
}
template<int M> void checkMod() {
	using Z = ModInt<M>;
	assert(Z().v == 0 && Z::mod() == M);
	vector<ll> edges{LLONG_MIN, INT_MIN, -2LL*M, -M, -1,
		0, 1, M-1, M, 2LL*M, INT_MAX, LLONG_MAX};
	for (ll a : edges) for (ll b : edges) checkPair<M>(a, b);
	if constexpr (M <= 17)
		for (int a = 0; a < M; ++a) for (int b = 0; b < M; ++b)
			checkPair<M>(a, b);
	for (ll a : edges) {
		for (auto e : {0ULL, 1ULL, 2ULL, 63ULL, ULLONG_MAX})
			assert(pow(Z(a), e).v == power<M>(reduce<M>(a), e));
		assert(pow(Z(a), LLONG_MAX).v == power<M>(reduce<M>(a), LLONG_MAX));
		int r = 1 % M;
		for (int e = 0; e < 40; ++e) {
			assert(pow(Z(a), e).v == r);
			r = reduce<M>(i128(r)*reduce<M>(a));
		}
	}
	mt19937_64 rng(123);
	for (int it = 0; it < 2000; ++it) {
		ll x = bit_cast<ll>(rng()), y = bit_cast<ll>(rng());
		checkPair<M>(x, y);
		Z a(x), b(y), c(rng());
		assert((a+b)+c == a+(b+c));
		assert((a*b)*c == a*(b*c));
		assert(a*(b+c) == a*b+a*c);
	}
	for (auto a : {0U, 1U, UINT_MAX}) assert(Z(a).v == int(a % M));
	for (auto a : {0ULL, 1ULL, ULLONG_MAX}) assert(Z(a).v == int(a % M));
	assert(Z(false) == 0 && Z(true) == 1);
	for (auto a : {numeric_limits<i128>::min(), numeric_limits<i128>::max()})
		assert(Z(a).v == reduce<M>(a));
	u128 big = numeric_limits<u128>::max();
	assert(Z(big).v == int(big % M));
	assert(pow(Z(2), big).v == power<M>(2 % M, big));
	for (int i = SHRT_MIN; i <= SHRT_MAX; ++i)
		assert(Z(static_cast<short>(i)).v == reduce<M>(i));
	for (int i = SCHAR_MIN; i <= SCHAR_MAX; ++i)
		assert(Z(static_cast<signed char>(i)).v == reduce<M>(i));
	for (int i = CHAR_MIN; i <= CHAR_MAX; ++i)
		assert(Z(static_cast<char>(i)).v == reduce<M>(i));
	assert(Z(static_cast<unsigned char>(UCHAR_MAX)).v == UCHAR_MAX % M);
	assert(Z(static_cast<unsigned short>(USHRT_MAX)).v == USHRT_MAX % M);
	assert(Z(LONG_MIN).v == reduce<M>(LONG_MIN));
	assert(Z(ULONG_MAX).v == int(ULONG_MAX % M));
	assert(Z::raw(M-1) == Z(M-1));
	stringstream input("-1 " + to_string(LLONG_MIN) + " " + to_string(LLONG_MAX));
	Z a;
	input >> a; assert(a == -1);
	input >> a; assert(a == LLONG_MIN);
	input >> a; assert(a == LLONG_MAX);
	Z saved = a;
	input >> a; assert(input.fail() && a == saved);
	stringstream invalid("invalid");
	invalid >> a; assert(invalid.fail() && a == saved);
	stringstream overflow("9223372036854775808");
	overflow >> a; assert(overflow.fail() && a == saved);
	ostringstream output;
	const Z negative(-1);
	output << negative;
	assert(output.str() == to_string(M-1));
}
void rejects(auto f) {
	pid_t child = fork();
	assert(child >= 0);
	if (!child) {
		close(STDERR_FILENO);
		alarm(2);
		f();
		_Exit(0);
	}
	int status;
	assert(waitpid(child, &status, 0) == child);
	assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
}
int main() {
	checkMod<1>(); checkMod<2>(); checkMod<12>(); checkMod<17>();
	checkMod<MOD>(); checkMod<1000000007>();
	checkMod<1224736769>(); checkMod<INT_MAX>();
	using Z = ModInt<1073741824>;
	assert(Z(1073741823) + Z(1073741823) == 1073741822);
	vector<mint> v{1, 2, 3};
	assert(find(all(v), mint(2)) == v.begin()+1);
	assert(v == vector<mint>({1, 2, 3}));
	vector<vector<mint>> matrix{{0, 1}, {2, 3}};
	assert(gaussian::getRow(matrix, 2, 0, 0) == 1);
	rejects([] { (void)pow(mint(2), -1LL); });
	rejects([] { (void)pow(mint(2), LLONG_MIN); });
	rejects([] { (void)invert(mint(0)); });
	rejects([] { (void)(mint(1)/mint(0)); });
	rejects([] { (void)invert(ModInt<12>(2)); });
	rejects([] { (void)mint::raw(-1); });
	rejects([] { (void)mint::raw(MOD); });
	cout << "Tests passed!" << endl;
}
