/**
 * Author: Jacob
 * Date: back in the day
 * License: CC0
 * Source: me
 * Description: examples for PBDS BBST, mergeable heaps and rope.
 * Time: O(\log N)
 * Status: tested
 */
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
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
void pbds() {	
	Tree<int> t, t2; t.insert(8);
	auto it = t.insert(10).first;
	assert(it == t.lower_bound(9));
	assert(t.order_of_key(10) == 1);
	assert(t.order_of_key(11) == 2);
	assert(*t.find_by_order(0) == 8);
	t.join(t2); // assuming T < T2 or T > T2, merge t2 into t

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
	v.insert(v.mutable_begin() + 2, cur);
    v.insert(v.mutable_begin(), cur); //to start (2 3 4 1 5 6)
    // v.insert(v.mutable_reference_at(0), cur); // to ONE AFTER start (1 2 3 4 5 6)
    // v.insert(v.mutable_begin() + 2, cur); // to TWO AFTER start (1 5 2 3 4 6)
}