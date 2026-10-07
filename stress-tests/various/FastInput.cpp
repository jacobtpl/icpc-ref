#include "../utilities/template.h"
#include <unistd.h>

#include "../../content/various/FastInput.h"

// Compile with -DBENCH for timings against scanf / cin.
// gc() keeps its buffer in function-local statics, so every test below must consume its
// input up to and including EOF (gc() returning 0) before stdin is reopened.
const int BUF_SIZE = 1 << 16;

string tempdirname;
string tempfilename;

void openInput(const string& s) {
	ofstream fout(tempfilename, ios::binary);
	fout << s;
	fout.close();
	FILE* ret = freopen(tempfilename.c_str(), "r", stdin);
	assert(ret == stdin);
}

// s must not contain 0 bytes
void testChars(const string& s) {
	openInput(s);
	for (char c : s) assert(gc() == c);
	rep(i,0,3) assert(gc() == 0);
}

void test(const string& s, vi ints) {
	openInput(s);
	for (int x : ints) {
		int y = readInt();
		if (x != y) {
			cerr << "On input " << s.substr(0, 100) << ", read " << y << " but expected " << x << endl;
		}
		assert(x == y);
	}
	while (gc() != 0);
	rep(i,0,3) assert(gc() == 0);
}

int main() {
	char pattern[] = "/tmp/fastinputXXXXXX";
	tempdirname = mkdtemp(pattern);
	tempfilename = tempdirname + "/stdin.txt";
	mt19937 rng(99);

	// First test that the getchar implementation is correct:
	testChars("");
	testChars("a");
	testChars("ab");
	string s;
	for (int i = 0; i < BUF_SIZE; i++) s += (char)(i % 13 + 1);
	testChars(s);
	for (int i = 0; i < BUF_SIZE * 10 + 1; i++) s += (char)(i % 255 + 1); // includes chars >= 128
	testChars(s);
	for (int i = 0; i < BUF_SIZE - 2; i++) s += (char)(i % 13 + 1);
	testChars(s);
	for (int len : {BUF_SIZE - 1, BUF_SIZE + 1, 2 * BUF_SIZE, 3 * BUF_SIZE - 1})
		testChars(string(len, 'x'));
	for (int i = 0; i < BUF_SIZE + 2; i++) {
		assert(gc() == 0);
	}

	// Then test that readInt() is:
	test("1", {1});
	test("12", {12});
	test("9\n", {9});
	test("12 ", {12});
	test("-23\n", {-23});
	test(" -4", {-4});
	test(" 5\n", {5});
	test("1 -2 ", {1, -2});
	test("  -34   56   ", {-34, 56});
	test(" \t\r\n5 -2 ", {5});
	test("0", {0});
	test("-0 0 007 -007\r\n", {0, 0, 7, -7});
	test("1\r\n2\r\n3\r\n", {1, 2, 3});
	test("1\t2\n\n\n3", {1, 2, 3});
	test("1000000007 -1000000007 999999999", {1000000007, -1000000007, 999999999});
	// Largest values (INT_MIN itself cannot be read: -readInt() would negate 2^31).
	test("2147483647", {INT_MAX});
	test("-2147483647 2147483647\n", {-INT_MAX, INT_MAX});
	test("2147483599 2147483600 2147483639 2147483640 2147483646 -2147483640",
		{2147483599, 2147483600, 2147483639, 2147483640, 2147483646, -2147483640});

	// Numbers straddling the buffer boundary, for every offset.
	rep(pad,0,30) {
		string t(BUF_SIZE - 15 + pad, ' ');
		t += "-1234567890 2147483647\n-7";
		test(t, {-1234567890, INT_MAX, -7});
	}

	// Random files: random magnitudes and random whitespace, compared against what was written.
	rep(it,0,300) {
		int n = (int)(rng() % (it < 250 ? 200 : 60000)) + 1;
		vi v(n);
		string t;
		rep(i,0,n) {
			int bits = (int)(rng() % 32);
			ll x = bits == 31 ? (ll)(rng() % 2147483648u) : (ll)(rng() & ((1u << bits) - 1));
			if (rng() % 2) x = -x;
			v[i] = (int)x;
			int ws = rng() % 8 ? 0 : (int)(rng() % 4);
			if (i == 0 && rng() % 2) ws++;
			rep(j,0,ws) t += " \n\t\r"[rng() % 4];
			t += to_string(v[i]);
			t += " \n"[rng() % 2];
		}
		if (rng() % 2) t.pop_back(); // no trailing whitespace: number terminated by EOF
		test(t, v);
	}

#ifdef BENCH
	{
		const int N = 5'000'000;
		vi v(N);
		string t;
		rep(i,0,N) v[i] = (int)(rng() % 2000000001u) - 1000000000, t += to_string(v[i]), t += '\n';
		openInput(t);
		auto t0 = chrono::steady_clock::now();
		ll s1 = 0; rep(i,0,N) s1 += readInt();
		auto t1 = chrono::steady_clock::now();
		while (gc() != 0);
		FILE* f = fopen(tempfilename.c_str(), "r");
		ll s2 = 0; int x;
		rep(i,0,N) { int r = fscanf(f, "%d", &x); assert(r == 1); s2 += x; }
		fclose(f);
		auto t2 = chrono::steady_clock::now();
		ios::sync_with_stdio(0);
		ifstream fin(tempfilename);
		ll s3 = 0; rep(i,0,N) fin >> x, s3 += x;
		auto t3 = chrono::steady_clock::now();
		assert(s1 == s2 && s2 == s3 && s1 == accumulate(all(v), 0LL));
		cerr << "5e6 ints (" << sz(t) / 1000000 << " MB): readInt " << chrono::duration<double>(t1 - t0).count()
			<< " s, scanf " << chrono::duration<double>(t2 - t1).count()
			<< " s, ifstream>> " << chrono::duration<double>(t3 - t2).count() << " s" << endl;
	}
#endif

	unlink(tempfilename.c_str());
	rmdir(tempdirname.c_str());
	cout << "Tests passed!" << endl;
}
