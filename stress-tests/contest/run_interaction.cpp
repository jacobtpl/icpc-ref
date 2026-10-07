#include "../utilities/template.h"

// Tests content/contest/run_interaction.cpp by running real grader/user
// processes (shell scripts) through it.
namespace fs = std::filesystem;
string dir;

void script(const string& name, const string& body) {
	string p = dir + "/" + name;
	ofstream(p) << "#!/bin/bash\ncd '" << dir << "'\n" << body;
	fs::permissions(p, fs::perms::owner_all);
}
string slurp(const string& name) {
	stringstream ss; ss << ifstream(dir + "/" + name).rdbuf();
	string s = ss.str();
	while (!s.empty() && isspace(s.back())) s.pop_back();
	return s;
}
// runs "./ri args" with a 10 s limit; returns exit code (124 = timeout)
int run(const string& args) {
	int r = system(("cd '" + dir + "' && rm -f *.log && timeout 10 ./ri " + args +
		" > stdout.log 2> stderr.log").c_str());
	assert(WIFEXITED(r));
	return WEXITSTATUS(r);
}
int fails = 0;
void expect(bool ok, const string& what) {
	if (ok) return;
	if (fails++ < 10) cout << "FAILED: " << what << endl;
}

int main() {
	string src = fs::absolute(fs::path(__FILE__).parent_path() /
		"../../content/contest/run_interaction.cpp").string();
	assert(fs::exists(src));
	char tmpl[] = "/tmp/ritestXXXXXX";
	assert(mkdtemp(tmpl)); dir = tmpl;
	assert(system(("g++ -std=c++17 -O2 -Wall -Wextra -Wconversion '" + src +
		"' -o " + dir + "/ri").c_str()) == 0);

	// grader: number guessing; logs its argument count and what it saw
	script("grader.sh", R"(echo $# > grader_argc.log
n=0
while read -t 5 x; do
	n=$((n+1))
	if [ "$x" -lt $SECRET ]; then echo "<"
	elif [ "$x" -gt $SECRET ]; then echo ">"
	else echo "="; echo "ok $n" > grader.log; exit 0; fi
done
echo "eof $n" > grader.log
)");
	script("user.sh", R"(echo $# > user_argc.log
lo=0; hi=1000000; n=0
while true; do
	mid=$(( (lo+hi)/2 )); echo $mid; n=$((n+1))
	read -t 5 r || { echo "eof" > user.log; exit 3; }
	case "$r" in
		"<") lo=$((mid+1));; ">") hi=$((mid-1));;
		"=") echo "ok $mid $n" > user.log; exit 0;;
	esac
done
)");
	// 1. normal interaction, many random secrets
	mt19937 rng(3);
	rep(it,0,30) {
		int secret = it == 0 ? 0 : it == 1 ? 1000000 : (int)(rng() % 1000001);
		setenv("SECRET", to_string(secret).c_str(), 1);
		int rc = run("./grader.sh ./user.sh");
		expect(rc == 0, "normal run exit code " + to_string(rc));
		string u = slurp("user.log"), g = slurp("grader.log");
		expect(u.rfind("ok " + to_string(secret) + " ", 0) == 0, "user result: " + u);
		expect(g.rfind("ok", 0) == 0, "grader result: " + g);
		// programs must be started without extra arguments
		expect(slurp("grader_argc.log") == "0",
			"grader started with " + slurp("grader_argc.log") + " extra args");
		expect(slurp("user_argc.log") == "0",
			"user started with " + slurp("user_argc.log") + " extra args");
		expect(slurp("stdout.log") == "", "unexpected output");
		if (fails) break;
	}
	// 2. large volume through both pipes (more than the 64K pipe buffer)
	script("g_big.sh", R"(s=0
for i in $(seq 1 30000); do echo $i; read -t 5 x || exit 1; s=$((s+x)); done
echo "$s" > grader.log
)");
	script("u_big.sh", R"(while true; do
	read -t 5 x; s=$?
	[ $s -gt 128 ] && { echo "timeout" > user.log; exit 4; }
	[ $s -ne 0 ] && break
	echo $((2*x))
done
echo "eof" > user.log
)");
	int rc = run("./g_big.sh ./u_big.sh");
	expect(rc == 0, "big run exit code " + to_string(rc));
	expect(slurp("grader.log") == "900030000", "big run sum " + slurp("grader.log"));
	expect(slurp("user.log") == "eof", "big run: user did not see EOF");

	// 3. grader exits early: the user must see EOF on stdin, not hang
	script("g_quit.sh", "echo hello\nexit 0\n");
	script("u_wait.sh", R"(while true; do
	read -t 3 x; s=$?
	[ $s -gt 128 ] && { echo "timeout" > user.log; exit 4; }
	[ $s -ne 0 ] && break
done
echo "eof" > user.log
)");
	rc = run("./g_quit.sh ./u_wait.sh");
	expect(rc == 0 && slurp("user.log") == "eof",
		"user does not get EOF after grader exits (rc " + to_string(rc) +
		", user: " + slurp("user.log") + ")");

	// 4. user exits early: the grader must see EOF on stdin
	script("g_wait.sh", R"(while true; do
	read -t 3 x; s=$?
	[ $s -gt 128 ] && { echo "timeout" > grader.log; exit 4; }
	[ $s -ne 0 ] && break
done
echo "eof" > grader.log
)");
	script("u_quit.sh", "echo 5\nexit 0\n");
	rc = run("./g_wait.sh ./u_quit.sh");
	string g;
	rep(i,0,50) { // the grader may outlive run_interaction
		if ((g = slurp("grader.log")) != "") break;
		usleep(100000);
	}
	expect(rc == 0 && g == "eof",
		"grader does not get EOF after user exits (grader: " + g + ")");

	// 5. wrong usage: message on the terminal and failure exit code
	for (string args : {"", "./grader.sh", "./grader.sh ./user.sh extra"}) {
		rc = run(args);
		expect(rc != 0 && rc != 124, "'./ri " + args + "' exit code " + to_string(rc));
		expect(slurp("stdout.log").find("Usage") != string::npos,
			"'./ri " + args + "' prints no usage message");
	}
	fs::remove_all(dir);
	if (fails) return 1;
	cout << "Tests passed!" << endl;
}
