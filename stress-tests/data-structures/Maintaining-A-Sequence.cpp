#include "../utilities/template.h"

// Maintaining-A-Sequence.cpp is a full solution (reads stdin, writes stdout) to
// NOI 2005 "Maintaining a Sequence": INSERT / DELETE / MAKE-SAME / REVERSE /
// GET-SUM / MAX-SUM on a sequence. Run its main() on in-memory streams and
// compare against a brute force on a plain vector.
// The file contains a second, index-based variant that is compiled out by its
// first line (#define USE_POINTER). To test that one as well:
//   sed 's/^#define USE_POINTER//' ../../content/data-structures/Maintaining-A-Sequence.cpp > /tmp/mas_index.cpp
//   g++ -std=c++17 -O2 '-DSOLUTION="/tmp/mas_index.cpp"' Maintaining-A-Sequence.cpp && ./a.out
static FILE *IN, *OUT;
#ifndef SOLUTION
#define SOLUTION "../../content/data-structures/Maintaining-A-Sequence.cpp"
#endif
namespace sol {
#define scanf(...) fscanf(IN, __VA_ARGS__)
#define printf(...) fprintf(OUT, __VA_ARGS__)
#define main sol_main
#include SOLUTION
#undef main
#undef scanf
#undef printf
}

mt19937 rng(12345);
int ri(int a, int b) { return uniform_int_distribution<int>(a, b)(rng); }

string run(string& in) {
	char* buf = 0; size_t len = 0;
	IN = fmemopen(in.data(), in.size(), "r");
	OUT = open_memstream(&buf, &len);
	sol::sol_main();
	fclose(IN); fclose(OUT);
	string res(buf, len);
	free(buf);
	return res;
}

// maxLen: cap on sequence length, V: |value| bound, allowEmpty: generate
// zero-length ranges (tot = 0) and let the sequence become empty.
void test(int n, int m, int maxLen, int V, bool allowEmpty) {
	int lo = ri(0, 3) ? -V : (ri(0, 1) ? -V : 1), hi = ri(0, 5) ? V : -1;
	if (lo > hi) lo = -V, hi = V;
	vi a(n);
	string in = to_string(n) + " " + to_string(m) + "\n", exp;
	for (int& x : a) x = ri(lo, hi), in += to_string(x) + " ";
	in += "\n";
	int mn = allowEmpty ? 0 : 1;
	rep(it,0,m) {
		int len = sz(a), op = ri(0, 5);
		if (op == 0 && len < maxLen) {
			int p = ri(0, len), t = ri(mn, min({maxLen - len, 2 * n + 2, 300}));
			vi v(t);
			in += "INSERT " + to_string(p) + " " + to_string(t);
			for (int& x : v) x = ri(lo, hi), in += " " + to_string(x);
			a.insert(a.begin() + p, all(v));
		} else if (op == 5 || len < 1 + !allowEmpty) {
			if (!len) { it--; continue; }
			in += "MAX-SUM";
			ll best = LLONG_MIN, cur = 0;
			for (int x : a) cur = max(cur, 0LL) + x, best = max(best, cur);
			exp += to_string(best) + "\n";
		} else {
			int t = ri(mn, len - (op == 1 && !allowEmpty)), p = ri(1, len - t + 1);
			if (!ri(0, 9)) t = len - (op == 1 && !allowEmpty), p = 1;
			auto b = a.begin() + (p - 1), e = b + t;
			string s = " " + to_string(p) + " " + to_string(t);
			if (op == 1) in += "DELETE" + s, a.erase(b, e);
			else if (op == 2) {
				int c = ri(lo, hi);
				in += "MAKE-SAME" + s + " " + to_string(c), fill(b, e, c);
			}
			else if (op == 3) in += "REVERSE" + s, reverse(b, e);
			else {
				in += "GET-SUM" + s;
				exp += to_string(accumulate(b, e, 0LL)) + "\n";
			}
		}
		in += "\n";
	}
	string got = run(in);
	if (got != exp) {
		cerr << "MISMATCH on input:\n" << in << "expected:\n" << exp << "got:\n" << got;
		abort();
	}
}

int main(int argc, char** argv) {
	if (argc > 1) { // benchmark mode: ./a.out bench
		for (int V : {1000, 1}) {
			auto t0 = chrono::steady_clock::now();
			test(200000, 200000, 500000, V, 0);
			cerr << "N=M=2e5, |v|<=" << V << ": " << chrono::duration<double>(
				chrono::steady_clock::now() - t0).count() << " s (incl. brute force)\n";
		}
		// Problem limits: 5e5 elements at a time, 2e4 ops, 4e6 inserted in total.
		string in = "500000 60000\n";
		rep(i,0,500000) in += to_string(ri(-1000, 1000)) + " ";
		in += "\n";
		rep(i,0,20000) {
			if (i % 2) {
				in += "INSERT " + to_string(ri(0, 499800)) + " 200";
				rep(j,0,200) in += " " + to_string(ri(-1000, 1000));
				in += "\nMAX-SUM\n";
			} else {
				in += "REVERSE " + to_string(ri(1, 250000)) + " 250000\n";
				in += "MAKE-SAME " + to_string(ri(1, 499000)) + " 1000 " + to_string(ri(-1000, 1000)) + "\n";
				in += "GET-SUM " + to_string(ri(1, 250000)) + " 250000\n";
				in += "DELETE " + to_string(ri(1, 499800)) + " 200\n";
			}
		}
		auto t0 = chrono::steady_clock::now();
		string out = run(in);
		cerr << "N=5e5, 6e4 ops, 2e6 inserted: " << chrono::duration<double>(
			chrono::steady_clock::now() - t0).count() << " s\n";
		return 0;
	}
	rep(it,0,20000) {
		bool emp = it % 2;
		test(ri(!emp, 6), ri(1, 30), 12, ri(0, 3) ? 3 : 1000, emp);
	}
	rep(it,0,1000) test(ri(1, 200), 300, 400, ri(0, 1) ? 2 : 1000, it % 2);
	rep(it,0,10) test(ri(1, 3000), 3000, 5000, 1000, it % 2);
	test(100000, 20000, 500000, 1000, 0);
	cout << "Tests passed!" << endl;
}
