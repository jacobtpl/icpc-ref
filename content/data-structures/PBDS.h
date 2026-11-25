/**
 * Author: Jacob
 * Date: back in the day
 * License: CC0
 * Source: me
 * Description: examples for PBDS BBST, mergeable heaps and rope.
 * Time: O(log N), except split before GCC 16: linear in the size of the part split off.
 * Status: tested
 */
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/trie_policy.hpp>
#include <ext/pb_ds/priority_queue.hpp>
#include <ext/rope>
using namespace std;
using namespace __gnu_pbds;
using namespace __gnu_cxx;
template<class T>
using Tree = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template<class T>
using Heap = __gnu_pbds::priority_queue<T, less<T>, pairing_heap_tag>;
//binary_heap_tag, pairing_heap_tag, binomial_heap_tag, rc_binomial_heap_tag, thin_heap_tag

// Custom subtree aggregate: here a map int -> ll with the sum of values.
// operator() recomputes a node from its children (called on insert/erase).
// Changing a mapped value in place (t[k] = v, it->second = v) does NOT
// recompute: erase and re-insert instead.
template<class CIt, class It, class Cmp, class Alloc>
struct SumUpdate {
	typedef ll metadata_type; // aggregate stored per node
	static ll get(CIt it, CIt end) { return it == end ? 0 : it.get_metadata(); }
	void operator()(It it, CIt end) {
		const_cast<ll&>(it.get_metadata()) = (*it)->second
			+ get(it.get_l_child(), end) + get(it.get_r_child(), end);
	}
	virtual CIt node_begin() const = 0;
	virtual CIt node_end() const = 0;
	ll prefSum(int key) { // sum of values with key < given
		ll s = 0;
		for (CIt it = node_begin(), end = node_end(); it != end;) {
			if ((*it)->first < key) {
				s += (*it)->second + get(it.get_l_child(), end);
				it = it.get_r_child();
			} else it = it.get_l_child();
		}
		return s;
	}
};
using SumTree = tree<int, ll, less<int>, rb_tree_tag, SumUpdate>;
using Trie = trie<string, null_type, trie_string_access_traits<>,
	pat_trie_tag, trie_prefix_search_node_update>;

void pbds() {	
	Tree<int> t, t2; t.insert(8);
	auto it = t.insert(10).first;
	assert(it == t.lower_bound(9));
	assert(t.order_of_key(10) == 1);
	assert(t.order_of_key(11) == 2);
	assert(*t.find_by_order(0) == 8);
	t.join(t2); // assuming T < T2 or T > T2, merge t2 into t
	// join is O(log N); throws join_error if the key ranges overlap

	// split: keys > 3 move to b (old contents of b are dropped).
	// GCC 16+: O(log N). GCC <= 15: linear in the size of b afterwards.
	Tree<int> a, b;
	for (int x : {1, 2, 3, 4, 5, 6}) a.insert(x);
	a.split(3, b);
	assert(sz(a) == 3 && *b.begin() == 4 && b.order_of_key(6) == 2);

	// multiset: pair with a unique id
	Tree<pii> ms; int id = 0;
	for (int x : {5, 3, 5, 7}) ms.insert({x, id++});
	assert(ms.order_of_key({5, 0}) == 1); // # elements < 5
	assert(ms.order_of_key({6, 0}) == 3); // # elements <= 5
	assert(ms.find_by_order(2)->first == 5);
	ms.erase(ms.lower_bound({5, 0})); // erase one copy of 5

	SumTree st;
	st[1] = 10; st[5] = 20; st[9] = 30; // new keys: fine
	assert(st.prefSum(6) == 30);
	st.erase(5); st.insert({5, 1}); // change a value
	assert(st.prefSum(6) == 11);

	Trie tr;
	for (string s : {"apple", "app", "apt", "bat"}) tr.insert(s);
	auto [lo, hi] = tr.prefix_range("app"); // sorted: app, apple
	assert(distance(lo, hi) == 2 && *lo == "app");

	Heap<int> pq1, pq2;
	pq1.push(1); pq2.push(5);
	pq1.join(pq2); // merge pq2 into pq1
	assert(pq1.top()==5);
	auto pq_it = pq1.push(3);
	assert(pq1.top()==5);
	pq1.modify(pq_it,7); // modify-key in O(log N)
	assert(pq1.top()==7);
	
	int n=3;
    rope<int> v(n, 0);
    for (int i=0; i<n; i++) v.mutable_reference_at(i) = i + 1; // (1 2 3)
	for (int i=0; i<n; i++) v.push_back(i + n + 1); // (1 2 3 4 5 6)
    int l=1, r=3;
    rope<int> cur = v.substr(l, r-l+1); // 2 3 4
    v.erase(l, r-l+1); // 1 5 6
    v.insert(v.mutable_begin(), cur); //to start (2 3 4 1 5 6)
    // v.insert(1, cur); // to ONE AFTER start (1 2 3 4 5 6)
    // v.insert(v.mutable_begin() + 2, cur); // to TWO AFTER start (1 5 2 3 4 6)

	// Persistent array: copying a rope is O(1) (shares nodes),
	// later writes to either copy are O(log N) and don't affect the other.
	rope<int> old = v;
	v.replace(0, 9); // or v.mutable_reference_at(0) = 9
	assert(old[0] == 2 && v[0] == 9);
}