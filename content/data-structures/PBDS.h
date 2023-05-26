#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/priority_queue.hpp>
#include <ext/rope>
using namespace std;
using namespace __gnu_pbds;
using namespace __gnu_cxx;
#define pb push_back
int main() {
	tree<int,null_type,less<int>,rb_tree_tag,tree_order_statistics_node_update> t;
	t.insert(1);
	t.insert(7);
	t.insert(2);
	t.insert(5);
	assert(*t.find_by_order(2)==5);
	assert(t.find_by_order(5)==t.end());
	assert(t.order_of_key(4)==2);
	assert(t.order_of_key(5)==2);
	assert(t.order_of_key(8)==4);

	__gnu_pbds::priority_queue<int,less<int>,pairing_heap_tag> pq1,pq2;
	//binary_heap_tag, pairing_heap_tag, binomial_heap_tag, rc_binomial_heap_tag, thin_heap_tag
	pq1.push(1);
	pq2.push(5);
	pq1.join(pq2);
	assert(pq1.top()==5);
	auto it=pq1.push(3);
	assert(pq1.top()==5);
	pq1.modify(it,7);
	assert(pq1.top()==7);
	
	int n=5;
    rope <int> v(n, 0);
    for(int i = 0; i < n; ++i) v.mutable_reference_at(i) = i + 1; //push_back
    int l=1, r=3;
    rope <int> cur = v.substr(l, r-l+1);
    v.erase(l, r-l+1);
    v.insert(v.mutable_begin(), cur); //to start
    v.insert(v.mutable_reference_at(0), cur); //to one after start
    for (int x:v) printf("%d\n", x);
}