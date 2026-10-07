// Treap-jacob.h is a complete program (POJ 3580 SuperMemo): it has its own
// main() reading stdin. The test renames that main, feeds it generated input
// through stdin and compares what it prints with a vector-based oracle.
#include "../utilities/template.h"
#include <unistd.h>

// "int main() {...}" becomes a declaration plus "void treap_jacob_run() {...}"
// (the body has no return statement, so it cannot stay a non-void function)
#define main treap_jacob_unused(); void treap_jacob_run
#include "../../content/data-structures/Treap-jacob.h"
#undef main
#undef mid
#undef M
#undef ii

mt19937 rng(2718);
int rnd(int a, int b) { return a + (int)(rng() % (unsigned)(b - a + 1)); }

string in;
vector<int> expected;

void gen(int n, int q, int maxv) {
	vector<int> a(n);
	in += to_string(n) + "\n";
	for (int& x : a) x = rnd(-maxv, maxv), in += to_string(x) + " ";
	in += "\n" + to_string(q) + "\n";
	rep(it,0,q) {
		int len = sz(a), op = rnd(0, 5);
		if (len == 0 || op == 0 || (op == 1 && len < 3)) {
			int p = rnd(0, len), c = rnd(-maxv, maxv);
			in += "INSERT " + to_string(p) + " " + to_string(c) + "\n";
			a.insert(a.begin() + p, c); // after the p-th element
			continue;
		}
		int x = rnd(1, len), y = rnd(1, len);
		if (x > y) swap(x, y);
		if (rnd(0, 9) == 0) x = 1, y = len;
		if (op == 1) {
			in += "DELETE " + to_string(x) + "\n";
			a.erase(a.begin() + x - 1);
		} else if (op == 2) {
			int c = rnd(-maxv, maxv);
			in += "ADD " + to_string(x) + " " + to_string(y) + " " + to_string(c) + "\n";
			rep(i,x-1,y) a[i] += c;
		} else if (op == 3) {
			in += "REVERSE " + to_string(x) + " " + to_string(y) + "\n";
			reverse(a.begin() + x - 1, a.begin() + y);
		} else if (op == 4) {
			int c = rnd(0, 3) ? rnd(0, 2 * len) : rnd(-1000000, 1000000);
			in += "REVOLVE " + to_string(x) + " " + to_string(y) + " " + to_string(c) + "\n";
			int l = y - x + 1, k = ((c % l) + l) % l; // rotate right by k
			rotate(a.begin() + x - 1, a.begin() + y - k, a.begin() + y);
		} else {
			in += "MIN " + to_string(x) + " " + to_string(y) + "\n";
			expected.push_back(*min_element(a.begin() + x - 1, a.begin() + y));
		}
	}
	// read back the whole final sequence
	rep(i,0,sz(a)) {
		in += "MIN " + to_string(i + 1) + " " + to_string(i + 1) + "\n";
		expected.push_back(a[i]);
	}
}

int main() {
	string cases;
	// each case is "n, values, q, q ops" followed by a block of MIN i i queries; the
	// program reads q from the input, so the read-back queries are counted in q
	auto add = [&](int n, int q, int maxv) {
		in.clear();
		gen(n, q, maxv);
		// patch q: count the lines after the third one
		size_t p1 = in.find('\n'), p2 = in.find('\n', p1 + 1), p3 = in.find('\n', p2 + 1);
		ll lines = count(in.begin() + (ll)p3 + 1, in.end(), '\n');
		cases += in.substr(0, p2 + 1) + to_string(lines) + "\n" + in.substr(p3 + 1);
	};
	rep(it,0,30000) add(rnd(0, 8), rnd(0, 30), it % 2 ? 3 : 1000);
	rep(it,0,300) add(rnd(1, 300), rnd(0, 600), 100000);
	add(100000, 20000, 10000); // largest n the program supports
	add(1, 0, 5);

	char inPath[] = "/tmp/treap-jacob-in-XXXXXX", outPath[] = "/tmp/treap-jacob-out-XXXXXX";
	int fi = mkstemp(inPath), fo = mkstemp(outPath);
	assert(fi >= 0 && fo >= 0);
	assert(write(fi, cases.data(), cases.size()) == (ssize_t)cases.size());
	close(fi), close(fo);
	fflush(stdout);
	int savedOut = dup(1);
	assert(freopen(inPath, "r", stdin) && freopen(outPath, "w", stdout));
	treap_jacob_run();
	{ // the debug printer: "(v, mv to v[left] v[right]) " per node, in order
		int r;
		t.init();
		t.alloc(r, -5);
		t.output(r);
	}
	fflush(stdout);
	dup2(savedOut, 1);
	close(savedOut);

	FILE* f = fopen(outPath, "r");
	assert(f);
	size_t idx = 0;
	int x;
	while (fscanf(f, "%d", &x) == 1) {
		if (idx >= expected.size() || expected[idx] != x) {
			cerr << "answer #" << idx << ": got " << x << " expected "
				<< (idx < expected.size() ? to_string(expected[idx]) : "nothing") << endl;
			abort();
		}
		idx++;
	}
	char rest[100] = {};
	assert(fgets(rest, sizeof rest, f));
	fclose(f);
	if (string(rest) != "(-5, -5 to 0 0) ") {
		cerr << "output() printed \"" << rest << "\"" << endl;
		abort();
	}
	unlink(inPath), unlink(outPath);
	if (idx != expected.size()) {
		cerr << "got " << idx << " answers, expected " << expected.size() << endl;
		abort();
	}
	cout<<"Tests passed!"<<endl;
}
