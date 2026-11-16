#!/usr/bin/env bash
# Checks the 128-bit cin/cout operators in content/contest/template.cpp with
# every C++ compiler found on this machine, under several flag sets.
# Run from inside the repo:  bash stress-tests/various/check-128-io.sh
# Paste the whole output back if anything says FAIL.
# Works from the repo (any directory in it) or as a downloaded copy run there.
root=""
for d in "$(git rev-parse --show-toplevel 2>/dev/null)" "$(cd "$(dirname "$0")/../.." && pwd)"; do
	if [ -n "$d" ] && [ -f "$d/stress-tests/various/Template128IO.py" ] &&
		grep -q u128 "$d/content/contest/template.cpp" 2>/dev/null; then root=$d; break; fi
done
if [ -z "$root" ]; then
	echo "Cannot find the 128-bit I/O files. cd into your icpc-ref checkout, run"
	echo "  git fetch origin && git checkout devin/1794861157-template-128-io && git pull"
	echo "and run this script again from there."
	exit 1
fi
cd "$root" || exit 1
echo "repo: $root ($(git rev-parse --short HEAD 2>/dev/null))"
command -v python3 >/dev/null || { echo "python3 is required"; exit 1; }
echo "python: $(python3 --version 2>&1)"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/shim/bits"
printf '#include <%s>\n' algorithm cctype cstdlib cstring iomanip iostream \
	limits sstream string utility vector > "$tmp/shim/bits/stdc++.h"

echo "system: $(uname -srm)"
compilers=() seen=""
for c in clang++ g++ c++ g++-{10,11,12,13,14,15,16} \
	/opt/homebrew/bin/g++-{11,12,13,14,15,16} /usr/local/bin/g++-{11,12,13,14,15,16}; do
	command -v "$c" >/dev/null || continue
	v=$("$c" --version 2>/dev/null | head -1)
	case "$seen" in *"|$v|"*) continue;; esac
	seen="$seen|$v|"; compilers+=("$c")
done
[ ${#compilers[@]} -gt 0 ] || { echo "no C++ compiler found"; exit 1; }

flagsets=("" "-O2" "-std=c++11" "-std=c++17 -O2" "-std=gnu++17 -O2"
	"-std=c++20 -O2 -Wall -Wextra -Wconversion" "-std=gnu++20 -O0" "-Ofast")
fail=0
for c in "${compilers[@]}"; do
	echo; echo "== $c: $("$c" --version | head -1)"
	echo "   128-bit macros:$("$c" -x c++ -dM -E - </dev/null |
		grep -E 'SIZEOF_(INT128|FLOAT128)|(LDBL|FLT128)_MANT_DIG' | sed 's/#define//' | tr '\n' ';')"
	# Ask the compiler directly, not a macro: is __float128 a real 113-bit type?
	has=no
	echo 'int main(){__float128 a=1,e=1;for(int i=0;i<100;i++)e/=2;return !(sizeof a==16&&a+e!=a);}' \
		>"$tmp/q.cpp"
	"$c" "$tmp/q.cpp" -o "$tmp/q" 2>/dev/null && "$tmp/q" && has=yes
	echo "   compiler has real __float128: $has"
	shim=""
	if ! echo '#include <bits/stdc++.h>' | "$c" -x c++ -fsyntax-only - 2>/dev/null; then
		shim="-I$tmp/shim"
		echo "   no <bits/stdc++.h> here: using a temporary stand-in for this test"
	fi
	for f in "${flagsets[@]}"; do
		if out=$(CXX="$c" python3 stress-tests/various/Template128IO.py $shim $f 2>"$tmp/err"); then
			case "$has:$out" in
				yes:integers*) fail=1; r="FAIL (template fell back to long double)";;
				no:Template128IO*) fail=1; r="FAIL (unexpected __float128 path)";;
				*) r=PASS;;
			esac
			printf '   %s [%s] %s\n' "$r" "${f:-no flags}" "$out"
		else
			fail=1; printf '   FAIL [%s]\n' "${f:-no flags}"
			{ grep -E 'error|Error|assert|Abort|Errno|ld:|Undefined' "$tmp/err" | head -8
				echo "...last lines:"; tail -5 "$tmp/err"; } | cut -c1-300 | sed 's/^/      /'
		fi
	done
done
echo
[ $fail = 0 ] && echo "ALL PASSED" || echo "SOME FAILED"
exit $fail
