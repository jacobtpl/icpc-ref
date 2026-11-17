/**
 * Author: Simon Lindholm
 * Date: 2016-07-24
 * License: CC0
 * Source: Russian page
 * Description: Pre-computation of modular inverses: after initInv(), inv[i] = 1/i for 1 <= i < LIM. Assumes LIM <= mod and that mod is a prime.
 * Status: Works
 */
#pragma once

const ll mod = 1000000007, LIM = 200000;
ll inv[LIM];
void initInv() {
	inv[1] = 1;
	rep(i,2,LIM) inv[i] = mod - (mod / i) * inv[mod % i] % mod;
}
