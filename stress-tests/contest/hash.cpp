#include "../utilities/template.h"

// Tests content/contest/hash.sh: the hash must ignore whitespace and
// comments, and change when any token changes.
namespace fs = std::filesystem;
string dir, script;

string hashOf(const string& code) {
	string f = dir + "/in.cpp";
	ofstream(f) << code;
	string cmd = "bash '" + script + "' < '" + f + "' 2>" + dir + "/err";
	FILE* p = popen(cmd.c_str(), "r");
	assert(p);
	char buf[256]; string out;
	while (fgets(buf, sizeof buf, p)) out += buf;
	assert(pclose(p) == 0);
	while (!out.empty() && isspace(out.back())) out.pop_back();
	assert(sz(out) == 6);
	for (char c : out) assert(isxdigit(c));
	return out;
}

int main() {
	script = fs::absolute(fs::path(__FILE__).parent_path() /
		"../../content/contest/hash.sh").string();
	assert(fs::exists(script));
	char tmpl[] = "/tmp/hashtestXXXXXX";
	assert(mkdtemp(tmpl)); dir = tmpl;
	mt19937 rng(5);

	vector<string> toks = {"int", "main", "(", ")", "{", "x", "+=", "1",
		";", "}", "\"a b\"", "'c'", "/", "*", "<", "vi", "0x1f", "#"};
	rep(it,0,300) {
		int n = (int)(rng() % 30);
		vector<string> t;
		rep(i,0,n) t.push_back(toks[rng() % (sz(toks) - 1)]); // no '#'
		// variant A: one space between tokens; B: random whitespace
		// and comments between tokens
		string a, b;
		for (auto& s : t) {
			a += s + " ";
			b += s;
			int w = (int)(rng() % 6);
			if (w == 0) b += " ";
			if (w == 1) b += "\n\t";
			if (w == 2) b += " /* note */ ";
			if (w == 3) b += " // line comment\n";
			if (w == 4) b += "\r\n  ";
			if (w == 5) b += " /* multi\n line * / */\t";
		}
		a += "\n"; b += "\n";
		string ha = hashOf(a);
		assert(ha == hashOf(b));
		if (n) { // changing one token changes the hash
			string c;
			int at = (int)(rng() % n);
			rep(i,0,n) c += (i == at ? string("y") : t[i]) + " ";
			assert(hashOf(c + "\n") != ha);
		}
	}
	// edge cases
	string emptyHash = hashOf("");
	assert(emptyHash == "d41d8c"); // md5 of the empty string
	assert(hashOf("\n\n  \t\n") == emptyHash);
	assert(hashOf("// only a comment\n/* and\nanother */\n") == emptyHash);
	// comment markers inside string literals are not comments
	assert(hashOf("char*s=\"/*x*/\";\n") != hashOf("char*s=\"\";\n"));
	assert(hashOf("char*s=\"//x\";\n") != hashOf("char*s=\"\n"));
	// preprocessor directives are kept, not expanded
	assert(hashOf("#define A 1\nint x=A;\n") != hashOf("int x=A;\n"));
	assert(hashOf("#define A 1\nint x=A;\n") != hashOf("#define A 1\nint x=1;\n"));
	assert(hashOf("#define A 1\nint x = A; // c\n") == hashOf("#define A 1\nint x=A;\n"));
	assert(hashOf("#include <bits/stdc++.h>\nint main(){}\n") ==
		hashOf("#include <bits/stdc++.h>\n\nint main() {\n}\n"));
	// no trailing newline
	assert(hashOf("int x;") == hashOf("int x;\n"));
	// large input (~1.2 MB) stays fast
	string big;
	rep(i,0,60000) big += "int v" + to_string(i) + " = 0; // c\n";
	auto t0 = chrono::steady_clock::now();
	string hb = hashOf(big);
	double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
	assert(el < 5);
	big[10] = 'z';
	assert(hashOf(big) != hb);
	fs::remove_all(dir);
	cout << "Tests passed!" << endl;
}
