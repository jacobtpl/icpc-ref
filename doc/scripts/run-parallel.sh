#!/usr/bin/env bash
# Compile and run stress tests in parallel, one job per core.
# Usage: CXX=g++ CXXFLAGS="-O2 ..." run-parallel.sh test1.cpp test2.cpp ...
# A test passes if it compiles, exits 0 and, when REQUIRE_PASS_LINE=1, prints "tests passed".
set -u
: "${CXX:=g++}" "${CXXFLAGS:=-O2}" "${REQUIRE_PASS_LINE:=0}"
export CXX CXXFLAGS REQUIRE_PASS_LINE
out=$(mktemp -d)
export out
trap 'rm -rf "$out"' EXIT
ulimit -c 0

run_one() {
	name=$(echo "$1" | tr '/' '_')
	start=$SECONDS
	if $CXX -std=c++20 $CXXFLAGS "$1" -o "$out/$name" >"$out/$name.log" 2>&1 &&
		"$out/$name" >>"$out/$name.log" 2>&1 &&
		{ [ "$REQUIRE_PASS_LINE" != 1 ] || grep -qi "tests passed" "$out/$name.log"; }; then
		echo "ok   $1 ($((SECONDS - start))s)"
	else
		echo "FAIL $1 ($((SECONDS - start))s)"
		touch "$out/$name.failed"
	fi
}
export -f run_one

printf '%s\0' "$@" | xargs -0 -n 1 -P "$(getconf _NPROCESSORS_ONLN)" bash -c 'run_one "$0"'

status=0
for src in "$@"; do
	name=$(echo "$src" | tr '/' '_')
	if [ -e "$out/$name.failed" ]; then
		status=1
		echo "::group::FAILED $src"
		cat "$out/$name.log"
		echo "::endgroup::"
	elif [ -s "$out/$name.log" ]; then
		echo "::group::$src"
		cat "$out/$name.log"
		echo "::endgroup::"
	fi
done
exit $status
