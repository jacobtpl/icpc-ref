/**
 * Author: Siyong
 * Date: April 6, 2024
 * Description: Computes inv, fact, ifact. Requires prime MOD and MV <= MOD.
 *   Call init before use; choose requires 0 <= k <= n < MV.
 * Status: Tested many many times
 * Time: O(N)
 */
#include "ModInt.h"

int const MV = 2e6 + 10;
static_assert(MV <= MOD);
mint inv[MV], fact[MV], ifact[MV];
void init() {
	inv[1] = mint(1);
	for(int i = 2; i < MV; ++i)
		inv[i] = mint(MOD - MOD/i) * inv[MOD % i];
	fact[0] = ifact[0] = mint(1);
	for(int i = 1; i < MV; ++i) {
		fact[i] = mint(i) * fact[i-1];
		ifact[i] = inv[i] * ifact[i-1];
	}
}
mint choose(int n, int k) {
	assert(0 <= k && k <= n && n < MV); // or return 0
	return fact[n] * ifact[n-k] * ifact[k];
}
