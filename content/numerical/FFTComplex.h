/**
 * Author: Siyong
 * Date: 2022
 * License: CC0
 * Source: idk
 * Description: FFT but with complex numbers
 * Time: O(N \log N) with $N = |A|+|B|$ ($\tilde 1s$ for $N=2^{22}$)
 * Status: somewhat tested
 */
vector<C> conv_complex(const vector<C>& a, const vector<C>& b) {
	if (a.empty() || b.empty()) return {};
	vector<C> res(sz(a) + sz(b) - 1);
	int L = 32 - __builtin_clz(sz(res)), n = 1 << L;
	vector<C> in1(n), in2(n), out(n);
	copy(all(a), begin(in1));
	copy(all(b), begin(in2));
	fft(in1);
	fft(in2);
	rep(i,0,n) out[i] = conj(in1[i] * in2[i]);
	fft(out);
	rep(i,0,sz(res)) res[i] = conj(out[i]) / C(n, 0);
	return res;
}
