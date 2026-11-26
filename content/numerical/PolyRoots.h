/**
 * Author: Per Austrin
 * Date: 2004-02-08
 * License: CC0
 * Description: Finds the real roots of a polynomial in [xmin, xmax], sorted.
 *  A root up to 1e-9 (relative) outside the range is returned as the endpoint.
 *  Only roots where the sign changes are guaranteed: a root of even
 *  multiplicity (x^2 = 0) may be missed, and a multiple root may be repeated.
 *  Returns {} for a constant polynomial, including the zero polynomial.
 * Usage: polyRoots({{2,-3,1}},-1e9,1e9) // solve x^2-3x+2 = 0
 * Time: O(n^2 log(1/eps))
 * Status: stress-tested against polynomials built from known roots
 */
#pragma once

#include "Polynomial.h"

vector<double> polyRoots(Poly p, double xmin, double xmax) {
	while (sz(p.a) && p.a.back() == 0) p.a.pop_back();
	if (sz(p.a) < 2) return {};
	vector<double> ret;
	// Slack so that rounding does not lose a root at xmin or xmax.
	double d = 1e-9 * (1 + max(abs(xmin), abs(xmax)));
	if (sz(p.a) == 2) {
		double x = -p.a[0] / p.a[1];
		if (xmin - d <= x && x <= xmax + d) ret.push_back(clamp(x, xmin, xmax));
		return ret;
	}
	Poly der = p;
	der.diff();
	auto dr = polyRoots(der, xmin, xmax);
	dr.insert(dr.begin(), xmin - d);
	dr.push_back(xmax + d);
	rep(i,0,sz(dr)-1) {
		double l = dr[i], h = dr[i+1], fl = p(l), fh = p(h);
		if (min(fl, fh) <= 0 && max(fl, fh) >= 0) {
			bool sign = fl > fh; // p is monotone on [l, h]
			rep(it,0,60) { // while (h - l > 1e-8)
				double m = (l + h) / 2, f = p(m);
				if ((f <= 0) ^ sign) l = m;
				else h = m;
			}
			ret.push_back(clamp((l + h) / 2, xmin, xmax));
		}
	}
	return ret;
}
