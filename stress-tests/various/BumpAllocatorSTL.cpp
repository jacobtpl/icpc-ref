#include "../utilities/template.h"

#include "../../content/various/BumpAllocatorSTL.h"

mt19937 rng(4321);
int rnd() { return (int)(rng() >> 1); }

// Compile with -DBENCH for timings (small<> vs std::allocator).
struct Wide { __int128 x, y; };
__attribute__((noinline)) void fill(vector<Wide, small<Wide>>& v, int k) { for (auto& a : v) a.x = k, a.y = -k; }

template<class T> void checkAlign() {
	small<T> al;
	size_t n = rnd() % 5 + 1u;
	T* p = al.allocate(n);
	assert((size_t)p % alignof(T) == 0);
	assert(buf <= (char*)p && (char*)(p + n) <= buf + sizeof buf);
	memset((void*)p, 0x5a, n * sizeof(T));
	al.deallocate(p, n);
}

int main() {
	assert((size_t)buf % 16 == 0);

	// Raw allocations: aligned, inside buf, disjoint.
	char* prev = buf + sizeof buf;
	rep(it,0,20000) {
		checkAlign<char>(); checkAlign<short>(); checkAlign<int>();
		checkAlign<double>(); checkAlign<long double>(); checkAlign<Wide>();
		checkAlign<pair<char, ll>>();
		small<char> al;
		size_t n = rnd() % 37u;
		char* p = al.allocate(n);
		assert(p + n <= prev);
		prev = p;
	}

	// The usage from the header: adjacency lists.
	rep(it,0,200) {
		int n = rnd() % 50 + 1, m = rnd() % 300;
		vector<vector<int, small<int>>> ed(n);
		vector<vi> ref(n);
		rep(i,0,m) {
			int a = rnd() % n, b = rnd() % n;
			ed[a].push_back(b); ref[a].push_back(b);
		}
		rep(i,0,n) assert(equal(all(ed[i]), all(ref[i])));

		// Copy / move / swap of containers using the allocator.
		vector<int, small<int>> a = ed[0], b, c(ed[n-1]);
		b = a; assert(b == ed[0]);
		b = move(a); assert(b == ed[0]);
		swap(b, c); assert(c == ed[0] && b == ed[n-1]);
		ed[0] = ed[n-1]; assert(ed[0] == b);
		ed.push_back(c); ed.insert(ed.begin(), b);
		assert(ed[0] == b && ed.back() == c);
		rep(i,0,n) assert(equal(all(ed[i+1]), all(ref[i ? i : n-1])));
		vector<vector<int, small<int>>> ed2 = ed; ed2 = ed; assert(ed2 == ed);
	}

	// Other containers, checked against the default-allocator versions.
	rep(it,0,200) {
		map<int, ll, less<int>, small<pair<const int, ll>>> m; map<int, ll> rm;
		set<ll, less<ll>, small<ll>> s; set<ll> rs;
		deque<int, small<int>> dq; deque<int> rdq;
		list<int, small<int>> l, l2; list<int> rl;
		unordered_map<int, int, hash<int>, equal_to<int>, small<pair<const int, int>>> um;
		basic_string<char, char_traits<char>, small<char>> str; string rstr;
		priority_queue<int, vector<int, small<int>>> pq; priority_queue<int> rpq;
		rep(i,0,300) {
			int x = rnd() % 100;
			m[x] += i; rm[x] += i; s.insert(x); rs.insert(x);
			if (rnd() % 2) dq.push_back(x), rdq.push_back(x); else dq.push_front(x), rdq.push_front(x);
			l.push_back(x); l2.push_front(x); rl.push_back(x); rl.push_front(x);
			um[x]++; str += (char)('a' + x % 26); rstr += (char)('a' + x % 26);
			pq.push(x); rpq.push(x);
			if (rnd() % 4 == 0) { assert(pq.top() == rpq.top()); pq.pop(); rpq.pop(); }
		}
		l.splice(l.begin(), l2); assert(l2.empty());
		assert(equal(all(m), all(rm)) && equal(all(s), all(rs)) && equal(all(dq), all(rdq)));
		assert(equal(all(l), all(rl)) && equal(all(str), all(rstr)) && sz(um) == sz(rs));
		auto m2 = m; m2 = m; auto s2 = move(s); swap(m, m2);
		assert(equal(all(m2), all(rm)) && equal(all(s2), all(rs)));
		vector<Wide, small<Wide>> w(rnd() % 50 + 1); fill(w, it);
		for (auto& a : w) assert(a.x == it && a.y == -it);
		vector<bool, small<bool>> vb(rnd() % 200, true); vb.push_back(false);
		assert(count(all(vb), true) == sz(vb) - 1);
	}

#ifdef BENCH
	{
		const int N = 1'000'000, M = 5'000'000;
		vector<pii> es(M);
		for (auto& e : es) e = {rnd() % N, rnd() % N};
		auto t0 = chrono::steady_clock::now();
		vector<vector<int, small<int>>> ed(N);
		for (auto [a, b] : es) ed[a].push_back(b), ed[b].push_back(a);
		auto t1 = chrono::steady_clock::now();
		vector<vi> ref(N);
		for (auto [a, b] : es) ref[a].push_back(b), ref[b].push_back(a);
		auto t2 = chrono::steady_clock::now();
		rep(i,0,N) assert(equal(all(ed[i]), all(ref[i])));
		cerr << "adjacency lists N=1e6 M=5e6, small<int>: " << chrono::duration<double>(t1 - t0).count()
			<< " s, std::allocator: " << chrono::duration<double>(t2 - t1).count() << " s, used "
			<< ((sizeof buf - buf_ind) >> 20) << " MB of buf" << endl;
	}
#endif
	cout<<"Tests passed!"<<endl;
}
