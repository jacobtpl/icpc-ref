#include "../utilities/template.h"

// Tests content/contest/sy_template.cpp: it must compile with and
// without -DLOCAL, programs must behave identically in both modes, and
// with -DLOCAL out-of-range indexing must throw.
namespace fs = std::filesystem;
string dir, tpl;

const char* pre = R"(
#define main sy_main
#include SYT
#undef main
#include <algorithm>
#include <stdexcept>
#include <string>
)";

vector<string> progs = {
// 0: basic use, checked indexing
R"(vector<int> f() { vector<int> r(3, 7); return r; }
int main() {
	vector<vector<int>> g(3, vector<int>(4)); g[2][3] = 1;
	vector<int> a(5), b = a, c(all(a)); const vector<int> d = f();
	b = d; c = std::move(b); a[4] = -1; std::sort(all(a));
	vector<ll> e; e.push_back(1LL << 40); e.emplace_back(2);
	printf("%d %d %d %d %lld\n", a[0], (int)c.size(), d[2], g[2][3], e[0] + e[1]);
	int x = 3; bool r1 = ckmax(x, 5LL); bool r2 = ckmin(x, 7); bool r3 = ckmin(x, -2);
	printf("%d %d %d %d\n", r1, r2, r3, x);
	int t = 0;
#ifdef LOCAL
	try { e[2] = 3; } catch (std::out_of_range&) { t++; }
	try { t += d[3]; } catch (std::out_of_range&) { t++; }
	try { vector<int> z; t += z[0]; } catch (std::out_of_range&) { t++; }
	try { t += g[3][0]; } catch (std::out_of_range&) { t++; }
	try { t += g[0][4]; } catch (std::out_of_range&) { t++; }
#else
	t = 5;
#endif
	printf("%d\n", t);
})",
// 1: brace initialisation
R"(int main() {
	vector<int> v{5, 2}, w = {1, 2, 3}, u{4};
	vector<std::string> s{"a", "b"};
	vector<vector<int>> g = {{1, 2}, {3}};
	printf("%d %d %d %d %d %d %d\n", (int)v.size(), v[0], (int)w.size(), w[2],
		(int)u.size(), (int)s.size(), (int)g[0].size());
})",
// 2: vector<bool>
R"(int main() {
	vector<bool> b(3); b[1] = true; const vector<bool> c(2, true);
	printf("%d %d %d\n", (int)b[0], (int)b[1], (int)c[1]);
#ifdef LOCAL
	try { b[3] = 1; puts("no throw"); return 1; } catch (std::out_of_range&) {}
#endif
})",
};

// returns program output, or "COMPILE ERROR"
string run(int id, bool local) {
	string src = dir + "/p.cpp", exe = dir + "/p", out = dir + "/out";
	ofstream(src) << pre << progs[id];
	string cmd = "g++ -std=c++17 -O2 -w " + string(local ? "-DLOCAL " : "") +
		"-DSYT='\"" + tpl + "\"' " + src + " -o " + exe + " 2>/dev/null";
	if (system(cmd.c_str()) != 0) return "COMPILE ERROR";
	assert(system((exe + " > " + out).c_str()) == 0);
	stringstream ss; ss << ifstream(out).rdbuf();
	return ss.str();
}

int main() {
	tpl = fs::absolute(fs::path(__FILE__).parent_path() /
		"../../content/contest/sy_template.cpp").string();
	assert(fs::exists(tpl));
	char tmpl[] = "/tmp/sytestXXXXXX";
	assert(mkdtemp(tmpl)); dir = tmpl;
	int fails = 0;
	for (string d : {"", "-DLOCAL "}) {
		string cmd = "g++ -std=c++17 -O2 -Wall -Wextra -Wconversion " + d +
			"'" + tpl + "' -o " + dir + "/t";
		if (system(cmd.c_str()) != 0) {
			cout << "sy_template.cpp does not compile with flags: " << d << endl;
			fails++;
		}
	}
	rep(id,0,sz(progs)) {
		string a = run(id, 0), b = run(id, 1);
		if (a == "COMPILE ERROR" || a != b) {
			cout << "program " << id << ": without LOCAL:\n" << a
				<< "\nwith LOCAL:\n" << b << endl;
			fails++;
		}
	}
	fs::remove_all(dir);
	if (fails) return 1;
	cout << "Tests passed!" << endl;
}
