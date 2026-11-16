#!/usr/bin/env bash
# Checks the 128-bit cin/cout operators in content/contest/template.cpp with
# every C++ compiler found on this machine, under several flag sets.
# Run from anywhere:  bash stress-tests/various/check-128-io.sh
# Paste the whole output back if anything says FAIL.
cd "$(dirname "$0")/../.." || exit 1
command -v python3 >/dev/null || { echo "python3 is required"; exit 1; }
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
		grep -E 'SIZEOF_(INT128|FLOAT128)|LDBL_MANT_DIG' | sed 's/#define//' | tr '\n' ';')"
	shim=""
	if ! echo '#include <bits/stdc++.h>' | "$c" -x c++ -fsyntax-only - 2>/dev/null; then
		shim="-I$tmp/shim"
		echo "   no <bits/stdc++.h> here: using a temporary stand-in for this test"
	fi
	for f in "${flagsets[@]}"; do
		if out=$(CXX="$c" python3 stress-tests/various/Template128IO.py $shim $f 2>"$tmp/err"); then
			printf '   PASS [%s] %s\n' "${f:-no flags}" "$out"
		else
			fail=1; printf '   FAIL [%s]\n' "${f:-no flags}"
			grep -E 'error|Error|assert|Abort' "$tmp/err" | head -8 | sed 's/^/      /'
		fi
	done
done
echo
[ $fail = 0 ] && echo "ALL PASSED" || echo "SOME FAILED"
exit $fail
