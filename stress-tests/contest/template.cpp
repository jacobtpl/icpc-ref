#include "../utilities/template.h"

// Tests content/contest/template.cpp: it must compile as strict C++17
// without diagnostics, and its macros/helpers must behave as expected.
namespace fs = std::filesystem;

const char* driver = R"(
#define main template_main
#include TEMPLATE
#undef main
int calls = 0;
int bound() { calls++; return 5; }
int main() {
	mt19937 rng(1);
	// ckmin / ckmax against obvious oracle
	rep(it,0,200000) {
		int a = (int)(rng() % 21) - 10, b = (int)(rng() % 21) - 10;
		int x = a; bool r = ckmax(x, b);
		assert(x == max(a, b) && r == (b > a));
		x = a; r = ckmin(x, b);
		assert(x == min(a, b) && r == (b < a));
	}
	{ // mixed and extreme types
		ll x = LLONG_MIN; assert(ckmax(x, INT_MAX) && x == INT_MAX);
		assert(ckmax(x, LLONG_MAX) && x == LLONG_MAX && !ckmax(x, x));
		assert(ckmin(x, LLONG_MIN) && x == LLONG_MIN && !ckmin(x, 0));
		double d = 1.5; assert(ckmax(d, 2) && d == 2.0 && !ckmin(d, 2.5));
		string s = "b"; assert(ckmax(s, "c") && s == "c" && ckmin(s, string("a")));
		pii p = mp(1, 2); assert(ckmax(p, pii(1, 3)) && p.second == 3);
		vi v = {3, 1}; assert(ckmin(v[0], v[1]) && v[0] == 1);
	}
	// rep: half-open, empty when a >= b, parenthesised arguments
	int c = 0; rep(i,0,10) c += i; assert(c == 45);
	c = 0; rep(i,5,5) c++; rep(i,7,3) c++; assert(c == 0);
	c = 0; rep(i,1?2:9,1?4:0) c += i; assert(c == 5);
	c = 0; rep(i,-3,0) c += i; assert(c == -6);
	c = 0; rep(i,0,bound()) c++; assert(c == 5);
	vi v = {3, 1, 2}; v.pb(0);
	sort(all(v)); assert(v == (vi{0, 1, 2, 3}) && sz(v) == 4);
	int arr[] = {2, 1}; sort(all(arr)); assert(arr[0] == 1);
	vi e; assert(sz(e) == 0 && sz(e) - 1 < 0); // sz is signed
	assert(sizeof(ll) == 8);
	cout << "ok" << endl;
}
)";

int main() {
	string tpl = fs::absolute(fs::path(__FILE__).parent_path() /
		"../../content/contest/template.cpp").string();
	assert(fs::exists(tpl));
	char tmpl[] = "/tmp/tpltestXXXXXX";
	assert(mkdtemp(tmpl)); string dir = tmpl;
	string flags = "g++ -std=c++17 -O2 -Wall -Wextra -Wconversion ";
	// the template itself is valid ISO C++17 and warning-free
	string cmd = flags + "-Werror -pedantic-errors '" + tpl + "' -o " + dir + "/tpl";
	if (system(cmd.c_str()) != 0) {
		cout << "template.cpp does not compile cleanly as C++17" << endl;
		return 1;
	}
	assert(system((dir + "/tpl < /dev/null").c_str()) == 0);
	ofstream(dir + "/d.cpp") << driver;
	cmd = flags + "-w -DTEMPLATE='\"" + tpl + "\"' " + dir + "/d.cpp -o " + dir + "/d";
	assert(system(cmd.c_str()) == 0);
	assert(system((dir + "/d > " + dir + "/out").c_str()) == 0);
	string out; ifstream(dir + "/out") >> out;
	assert(out == "ok");
	fs::remove_all(dir);
	cout << "Tests passed!" << endl;
}
