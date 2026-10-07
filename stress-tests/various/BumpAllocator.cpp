#include "../utilities/template.h"

#include "../../content/various/BumpAllocator.h"

mt19937 rng(1234);
int rnd() { return (int)(rng() >> 1); }

// Compile with -DBENCH for timings (new vs malloc).
struct Node { Node *l, *r; int v; };
struct Wide { __int128 x, y; };
struct LD { long double x; };

// noinline so that the compiler has to trust the alignment of the pointers
__attribute__((noinline)) void fill(vector<Wide>& v, int k) { for (auto& a : v) a.x = k, a.y = -k; }
__attribute__((noinline)) void cp(Wide* a, Wide* b) { *a = *b; }

template<class T> T* chk(T* p) {
	assert((size_t)p % alignof(T) == 0);
	return p;
}

int main() {
	char* lo = buf; char* hi = buf + sizeof buf;

	// Random sized raw allocations: inside buf, pairwise disjoint, contents preserved.
	vector<pair<char*, size_t>> blocks;
	rep(it,0,20000) {
		size_t s = rnd() % 4 == 0 ? rnd() % 3000 : rnd() % 40u;
		char* p = (char*)operator new(s);
		assert(lo <= p && p + s <= hi);
		// operator new must return memory aligned for any object that fits in it
		for (size_t al = 16; al > 1; al /= 2) if (s >= al) { assert((size_t)p % al == 0); break; }
		rep(j,0,(int)s) p[j] = (char)(it * 31 + j);
		blocks.emplace_back(p, s);
	}
	rep(it,0,sz(blocks)) {
		auto [p, s] = blocks[it];
		rep(j,0,(int)s) assert(p[j] == (char)(it * 31 + j));
		if (it) assert(p + s <= blocks[it-1].first); // disjoint (bump goes downwards)
	}

	// Typed allocations interleaved with odd-sized ones.
	rep(it,0,20000) {
		string s(rnd() % 50u, 'a'); // odd sized allocation when > 15 chars
		chk(new char)[0] = 'x';
		*chk(new double) = 1.5;
		*chk(new ll) = it;
		chk(new LD)->x = 1;
		Node* n = chk(new Node{0, 0, it}); assert(n->v == it);
		Wide* a = chk(new Wide()); Wide* b = chk(new Wide()); b->x = it; cp(a, b);
		assert(a->x == it);
		int* arr = chk(new int[rnd() % 5 + 1]); arr[0] = 1;
		delete n; delete[] arr; // no-ops
	}

	// STL containers on top of the replaced operator new.
	rep(it,0,300) {
		string s(17u + rnd() % 3, 'a');
		vector<Wide> v(rnd() % 100 + 1); fill(v, it);
		for (auto& a : v) assert(a.x == it && a.y == -it);
		vector<__int128> w(100); iota(all(w), 0);
		__int128 t = 0; for (auto x : w) t += x * x;
		assert(t == 328350);
		map<int, ll> m; set<int> st; vi q;
		rep(i,0,200) { int x = rnd() % 1000; m[x] += i; st.insert(x); q.push_back(x); }
		sort(all(q)); q.erase(unique(all(q)), q.end());
		assert(sz(q) == sz(st) && sz(m) == sz(st) && equal(all(q), st.begin()));
	}

	// Zero-sized and large allocations.
	assert(operator new(0) != nullptr);
	char* big = (char*)operator new(100 << 20);
	assert(lo <= big && big + (100 << 20) <= hi);
	big[0] = big[(100 << 20) - 1] = 1;

#ifdef BENCH
	{
		const int N = 10'000'000;
		auto t0 = chrono::steady_clock::now();
		Node* head = 0;
		rep(i,0,N) head = new Node{head, 0, i};
		auto t1 = chrono::steady_clock::now();
		Node* h2 = 0;
		rep(i,0,N) { Node* x = (Node*)malloc(sizeof(Node)); *x = Node{h2, 0, i}; h2 = x; }
		auto t2 = chrono::steady_clock::now();
		ll s = 0; for (Node* x = head; x; x = x->l) s += x->v;
		for (Node* x = h2; x; x = x->l) s -= x->v;
		assert(s == 0);
		cerr << "1e7 x new Node (bump): " << chrono::duration<double>(t1 - t0).count() << " s, malloc: "
			<< chrono::duration<double>(t2 - t1).count() << " s" << endl;
	}
#endif
	cout<<"Tests passed!"<<endl;
}
