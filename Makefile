help:
	@echo "make test-compiles  - compile every header in content/ after the template"
	@echo "make test           - build and run every stress test in stress-tests/"
	@echo "Set CXX to pick the compiler (default g++), e.g. CXX=g++-16 CXXFLAGS='-include cassert' on macOS."

test:
	./doc/scripts/run-all.sh .

test-compiles:
	./doc/scripts/compile-all.sh .

.PHONY: help test test-compiles
