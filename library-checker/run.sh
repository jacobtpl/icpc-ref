#!/usr/bin/env bash
# Runs graph/<problem>_<impl>.cpp against the public Library Checker test data.
# Usage: ./run.sh <path to library-checker-problems> [solution.cpp ...]
# The test data must exist already: python3 generate.py -p <problem> in that repo.
set -u
cd "$(dirname "$0")"
LCP=$(realpath "$1"); shift
[ $# -eq 0 ] && set -- graph/*.cpp string/*.cpp
tmp=$(mktemp -d)
for sol in "$@"; do
	name=$(basename "$sol" .cpp)
	dir=
	for d in "$LCP"/*/*/; do
		p=$(basename "$d")
		if [[ $name == "$p" || $name == "$p"_* ]] && [ ${#d} -gt ${#dir} ]; then dir=$d; fi
	done
	if [ -z "$dir" ] || [ ! -x "$dir/checker" ]; then echo "$name: no test data, skipped"; continue; fi
	python3 bundle.py "$sol" > "$tmp/a.cpp"
	${CXX:-g++} ${CXXFLAGS:--O2 -std=gnu++20} "$tmp/a.cpp" -o "$tmp/a" || { echo "$name: compile error"; continue; }
	bad=0 worst=0 n=0
	for in in "$dir"/in/*.in; do
		c=$(basename "$in" .in)
		s=$(date +%s%N)
		timeout 20 "$tmp/a" < "$in" > "$tmp/out"; rc=$?
		ms=$(( ($(date +%s%N) - s) / 1000000 ))
		[ $ms -gt $worst ] && worst=$ms
		n=$((n + 1))
		if [ $rc -ne 0 ] || ! "$dir/checker" "$in" "$tmp/out" "$dir/out/$c.out" > /dev/null 2>&1; then
			bad=$((bad + 1)); echo "  $name: FAILED $c (exit $rc, $ms ms)"
		fi
	done
	echo "$name: $((n - bad))/$n passed, slowest $worst ms"
done
rm -rf "$tmp"
