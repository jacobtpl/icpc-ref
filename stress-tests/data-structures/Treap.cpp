#include "../utilities/template.h"

#include "../../content/data-structures/Treap.h"

pair<Node*, Node*> split2(Node* n, int v) {
	if (!n) return {};
	if (n->val >= v) {
		auto pa = split2(n->l, v);
		n->l = pa.second;
		n->recalc();
		return {pa.first, n};
	} else {
		auto pa = split2(n->r, v);
		n->r = pa.first;
		n->recalc();
		return {n, pa.second};
	}
}

// the example application from the header (commented out there)
void move(Node*& t, int l, int r, int k) {
	Node *a, *b, *c;
	tie(a,b) = split(t, l); tie(b,c) = split(b, r - l);
	if (k <= l) t = merge(ins(a, b, k), c);
	else t = merge(a, ins(c, b, k - r));
}

mt19937 rng(1234);
int rnd(int a, int b) { return a + (int)(rng() % (unsigned)(b - a + 1)); }

vi contents(Node* t) {
	vi v;
	each(t, [&](int x) { v.push_back(x); });
	return v;
}

// sequence with insert / remove / move / range add / range reverse / range min
void runSeq(int ops, int maxv) {
	deque<Node> nodes; // stable addresses
	Node* t = 0;
	vi a;
	rep(it,0,ops) {
		int op = rnd(0, 6), n = sz(a);
		if (op == 0 || (op <= 2 && n < 4)) {
			int pos = rnd(0, n), v = rnd(-maxv, maxv);
			nodes.emplace_back(v);
			t = ins(t, &nodes.back(), pos);
			a.insert(a.begin() + pos, v);
		} else if (op == 1) {
			if (!n) continue;
			int pos = rnd(0, n - 1);
			t = remove(t, pos);
			a.erase(a.begin() + pos);
		} else if (op == 2) {
			int i = rnd(0, n), j = rnd(0, n);
			if (i > j) swap(i, j);
			int k = rnd(0, n);
			if (i < k && k < j) continue;
			move(t, i, j, k);
			int nk = (k >= j ? k - (j - i) : k);
			vi iv(a.begin() + i, a.begin() + j);
			a.erase(a.begin() + i, a.begin() + j);
			a.insert(a.begin() + nk, all(iv));
		} else {
			int l = rnd(0, n), r = rnd(0, n);
			if (l > r) swap(l, r);
			Node *x, *y, *z;
			tie(x,y) = split(t, l); tie(y,z) = split(y, r - l); // y = [l, r)
			assert(cnt(x) == l && cnt(y) == r - l && cnt(z) == n - r);
			if (y) {
				if (op == 3) {
					int d = rnd(-maxv, maxv);
					ladd(y, d);
					rep(i,l,r) a[i] += d;
				} else if (op == 4) {
					y->rev ^= 1, swap(y->l, y->r);
					reverse(a.begin() + l, a.begin() + r);
				} else {
					int exp = *min_element(a.begin() + l, a.begin() + r);
					if (y->minval != exp) {
						cerr << "min of [" << l << "," << r << ") got " << y->minval
							<< " expected " << exp << endl;
						abort();
					}
				}
			}
			t = merge(merge(x, y), z);
		}
		assert(cnt(t) == sz(a));
		if (it % 5 == 0 || it == ops - 1) {
			vi v = contents(t);
			if (v != a) {
				cerr << "sequence mismatch after " << it + 1 << " ops" << endl;
				abort();
			}
		}
	}
}

int ra() {
	static unsigned x;
	x *= 4176481;
	x += 193861934;
	return x >> 1;
}

int main() {
	srand(3);
	rep(it,0,1000) {
		vector<Node> nodes;
		vi exp;
		rep(i,0,10) {
			nodes.emplace_back(i*2+2);
			exp.emplace_back(i*2+2);
		}
		Node* n = 0;
		rep(i,0,10)
			n = merge(n, &nodes[i]);

		int v = rand() % 25;
		int left = cnt(split2(n, v).first);
		int rleft = (int)(lower_bound(all(exp), v) - exp.begin());
		assert(left == rleft);
	}

	rep(it,0,10000) {
		vector<Node> nodes;
		vi exp;
		rep(i,0,10) nodes.emplace_back(i);
		rep(i,0,10) exp.emplace_back(i);
		Node* n = 0;
		rep(i,0,10)
			n = merge(n, &nodes[i]);

		int i = ra() % 11, j = ra() % 11;
		if (i > j) swap(i, j);
		int k = ra() % 11;
		if (i < k && k < j) continue;

		move(n, i, j, k);
		// cerr << i << ' ' << j << ' ' << k << endl;

		int nk = (k >= j ? k - (j - i) : k);
		vi iv(exp.begin() + i, exp.begin() + j);
		exp.erase(exp.begin() + i, exp.begin() + j);
		exp.insert(exp.begin() + nk, all(iv));

		int ind = 0;
		each(n, [&](int x) {
			// cerr << x << ' ';
			assert(x == exp[ind++]);
		});
		// cerr << endl;
	}

	{ // empty / single element
		assert(!split(0, 0).first && !split(0, 0).second && !merge(0, 0));
		Node a(5);
		assert(merge(&a, 0) == &a && merge(0, &a) == &a);
		assert(split(&a, 0).second == &a && split(&a, 1).first == &a);
		assert(remove(&a, 0) == 0 && contents(0).empty());
	}
	rep(it,0,30000) runSeq(rnd(1, 40), it % 2 ? 3 : 1000000);
	rep(it,0,200) runSeq(2000, 100000);
	cout<<"Tests passed!"<<endl;
}
