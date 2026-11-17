#!/usr/bin/env bash
# Builds and runs every stress test. CXX picks the compiler (default g++).
DIR=${1:-.}
CXX=${CXX:-g++}
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

declare -i pass=0
declare -i fail=0
failTests=""
ulimit -s 524288 2>/dev/null || ulimit -s hard 2>/dev/null # For 2-sat test
for test in $(find $DIR/stress-tests -name '*.cpp' | sort); do
    echo "$(basename $test): "
    SECONDS=0
    $CXX -Wall -Wfatal-errors -Wconversion -std=c++20 -O2 $CXXFLAGS $test -o $tmp/a.out && $tmp/a.out
    retCode=$?
    if (($retCode != 0)); then
        echo "Failed with $retCode"
        fail+=1
        failTests="$failTests$test\n"
    else
        pass+=1
    fi
    echo "Took $SECONDS seconds"
    echo
done
echo "$pass/$(($pass+$fail)) tests passed"
if (($fail == 0)); then
    echo "No tests failed"
    exit 0
else
    echo -e "These tests failed: \n$failTests"
    exit 1
fi
