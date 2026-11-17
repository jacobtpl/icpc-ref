#!/usr/bin/env bash
# Compiles every header in content/ on its own, after the contest template.
# CXX picks the compiler (default g++), CXXFLAGS adds flags.
DIR=${1:-.}
CXX=${CXX:-g++}
SCRIPT_DIR=$DIR/doc/scripts
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

echo "skipped: "
find $DIR/content -name '*.h' | grep -Ff $SCRIPT_DIR/skip_headers
declare -i pass=0
declare -i fail=0
failHeaders=""
for test in $(find $DIR/content -name '*.h' | grep -vFf $SCRIPT_DIR/skip_headers | sort); do
    if [[ $($CXX -dumpmachine) != x86_64* ]] && grep -q "x86-64 only" $test; then
        echo "skipped (x86-64 only): $test"; continue
    fi
    printf '#include "%s"\n#include "%s"\n' \
        "$(cd $DIR && pwd)/content/contest/template.cpp" "$(cd "$(dirname $test)" && pwd)/$(basename $test)" > $tmp/t.cpp
    if $CXX -Wall -Wextra -Wfatal-errors -Wconversion -std=c++20 $CXXFLAGS -fsyntax-only $tmp/t.cpp; then
        pass+=1
    else
        echo "FAILED: $test"
        fail+=1
        failHeaders="$failHeaders$test\n"
    fi
done
echo "$pass/$(($pass+$fail)) headers compile"
if (($pass == 0)); then
    echo "No headers found (make sure skip_headers doesn't have whitespace lines)"
    exit 1
elif (($fail == 0)); then
    exit 0
else
    echo -e "These headers failed: \n$failHeaders"
    exit 1
fi
