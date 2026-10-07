#include "../utilities/template.h"

#include "../../content/various/SmallPtr.h"

struct Node {
	int val;
	ptr<Node> next, l, r;
	Node(int v = 0) : val(v) {}
};

mt19937 rng(99);
int rnd(int lo, int hi) { return uniform_int_distribution<int>(lo, hi)(rng); }

ptr<Node> insert(ptr<Node> t, int v) {
	if (!t) return new Node(v);
	if (v < t->val) t->l = insert(t->l, v);
	else t->r = insert(t->r, v);
	return t;
}
void walk(ptr<Node> t, vi& out) {
	if (!t) return;
	walk(t->l, out); out.push_back(t->val); walk(t->r, out);
}

int main() {
	assert(sizeof(ptr<Node>) == 4);
	assert(sizeof(Node) == 16);
	// null
	ptr<Node> nul;
	assert(!nul);
	assert(!ptr<Node>(nullptr));
	// round trips against raw pointers
	vector<Node*> raw;
	vector<ptr<Node>> sm;
	rep(i,0,200000) {
		Node* p = new Node(i);
		raw.push_back(p);
		sm.emplace_back(p);
	}
	rep(i,0,sz(raw)) {
		assert(sm[i]);
		assert(&*sm[i] == raw[i]);
		assert(sm[i]->val == i && (*sm[i]).val == i);
		assert(&sm[i][0] == raw[i]);
	}
	// arrays and operator[]
	rep(it,0,2000) {
		int n = rnd(1, 50);
		ll* arr = new ll[n];
		ptr<ll> p = arr;
		rep(i,0,n) p[i] = (ll)i * i - it;
		rep(i,0,n) assert(arr[i] == (ll)i * i - it && &p[i] == arr + i);
	}
	// chars: every byte offset must be representable
	rep(it,0,2000) {
		char* c = new char[rnd(1, 7)];
		ptr<char> p = c;
		assert(p && &*p == c);
	}
	// linked list
	{
		ptr<Node> head;
		rep(i,0,300000) {
			ptr<Node> x = new Node(i);
			x->next = head; head = x;
		}
		int e = 300000;
		for (ptr<Node> x = head; x; x = x->next) assert(x->val == --e);
		assert(e == 0);
	}
	// BST against std::multiset
	rep(it,0,300) {
		ptr<Node> root; multiset<int> ms;
		int n = rnd(0, 300);
		rep(i,0,n) { int v = rnd(-50, 50); root = insert(root, v); ms.insert(v); }
		vi out; walk(root, out);
		assert(out == vi(all(ms)));
	}
	// a big allocation that pushes the bump pointer far down
	{
		char* big = new char[400 << 20];
		ptr<char> p = big;
		assert(p && &*p == big);
		big[0] = 1; big[(400 << 20) - 1] = 2;
		assert(p[0] == 1 && p[(400 << 20) - 1] == 2);
		Node* q = new Node(5);
		ptr<Node> pq = q;
		assert(pq && pq->val == 5 && &*pq == q);
	}
	cout << "Tests passed!" << endl;
}
