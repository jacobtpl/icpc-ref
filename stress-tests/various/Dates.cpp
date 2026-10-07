#include "../utilities/template.h"

#include "../../content/various/Dates.h"

// Compile with -DBENCH for timings.
bool leap(int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; }
int daysIn(int m, int y) {
	const int d[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	return d[m-1] + (m == 2 && leap(y));
}

int main() {
	// The example from the header.
	{
		int jd = dateToInt(3, 24, 2004), m, d, y;
		intToDate(jd, m, d, y);
		assert(jd == 2453089 && m == 3 && d == 24 && y == 2004 && intToDay(jd) == "Wed");
	}
	// Known anchors (Julian day number at noon of that date).
	assert(dateToInt(1, 1, 2000) == 2451545 && intToDay(2451545) == "Sat");
	assert(dateToInt(1, 1, 1970) == 2440588 && intToDay(2440588) == "Thu");
	assert(dateToInt(11, 17, 1858) == 2400001 && intToDay(2400001) == "Wed");
	assert(dateToInt(10, 15, 1582) == 2299161 && intToDay(2299161) == "Fri");
	assert(dateToInt(11, 24, -4713) == 0 && intToDay(0) == "Mon"); // JD 0 (proleptic Gregorian)
	assert(dateToInt(2, 29, 2000) + 1 == dateToInt(3, 1, 2000));
	assert(dateToInt(2, 28, 1900) + 1 == dateToInt(3, 1, 1900));
	assert(dateToInt(12, 31, 9999) == 5373484);

	// Oracle: walk the (proleptic) Gregorian calendar one day at a time, starting from JD 0,
	// far past the documented 4-digit years. Every date must map to consecutive integers,
	// the inverse must round-trip and the weekday must advance cyclically.
	int y = -4713, m = 11, d = 24, jd = 0, dow = 0;
	ll cnt = 0;
	while (y <= 300000) {
		assert(dateToInt(m, d, y) == jd);
		int m2, d2, y2;
		intToDate(jd, m2, d2, y2);
		assert(m2 == m && d2 == d && y2 == y);
		assert(intToDay(jd) == dayOfWeek[dow]);
		if (++d > daysIn(m, y)) {
			d = 1;
			if (++m == 13) m = 1, y++;
		}
		jd++; dow = (dow + 1) % 7; cnt++;
	}
	assert(cnt > 110'000'000);

	// Independent oracle: the C library, for every day of 1900..2400.
	for (ll t = -2208988800LL; t < 13569465600LL; t += 86400) {
		time_t tt = (time_t)t;
		tm* g = gmtime(&tt);
		int j = dateToInt(g->tm_mon + 1, g->tm_mday, g->tm_year + 1900);
		assert(j == 2440588 + t / 86400);
		assert(intToDay(j) == dayOfWeek[(g->tm_wday + 6) % 7]);
	}

#ifdef BENCH
	{
		auto t0 = chrono::steady_clock::now();
		ll s = 0;
		rep(it,0,10) rep(j,1721426,5373485) { // 1/1/1 .. 12/31/9999, ten times
			int a, b, c;
			intToDate(j, a, b, c);
			s += dateToInt(a, b, c) - j;
		}
		auto t1 = chrono::steady_clock::now();
		assert(s == 0);
		cerr << "3.65e7 x (intToDate + dateToInt): " << chrono::duration<double>(t1 - t0).count() << " s" << endl;
	}
#endif
	cout<<"Tests passed!"<<endl;
}
