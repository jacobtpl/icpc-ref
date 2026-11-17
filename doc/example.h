/**
 * Author: Chen Xing
 * Date: 2009-03-13
 * License: CC0
 * Source: Conversation with Håkan and Ulf.
 * Description: Example of the file layout. Comments are plain text:
 *  write 2^k, a <= b, O(N log N).
 * Usage:
 *  int res = exampleFunction1(mice, cats, catsCount);
 * Time: O(M * N^4) where M = cats, N = mice
 * Status: Untested
 */
#pragma once

// Macros and typedefs from contest/template.cpp are available.
// Paste any reference file included here before this one.
#include "Example.h"

struct ExampleStructure {
	ll alpha;
	pii nodes;
	vi edges;
};

int exampleFunction1(vi mice, int* cats, int catsCount) {
	int sum = 0;
	rep(i, 1, catsCount)
		for (auto it : mice) sum += cats[i] * it;
	return sum;
}
