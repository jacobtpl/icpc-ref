/**
 * Author: someone on Codeforces
 * Date: 2017-03-14
 * Source: folklore
 * Description: A short self-balancing tree. It acts as a
 *  sequential container with log-time splits/joins, and
 *  is easy to augment with additional data. 0-index.
 * Time: $O(\log N)$
 * Status: stress-tested
 */
#pragma once

struct Node {
	Node *l = 0, *r = 0;
	int val, y, c = 1;
	Node(int val) : val(val), y(rand()) {}
	void recalc();
	void pushdown();
};
int cnt(Node* n) { return n ? n->c : 0; }
void ladd(Node *n, int add) {n->val+=add;n->minval+=add;n->lazyadd+=add;}
void Node::recalc() { c = cnt(l) + cnt(r) + 1;} // update range vals if needed
void Node::pushdown() {
	if (lazyadd) {
		if (l) ladd(l, lazyadd);
		if (r) ladd(r, lazyadd);
		lazyadd=0;}
	if (rev) {
		if (l) l->rev ^= 1, swap(l->l, l->r);
		if (r) r->rev ^= 1, swap(r->l, r->r);
		rev=0;}
}
template<class F> void each(Node* n, F f) {
	if (n) { each(n->l, f); f(n->val); each(n->r, f); }
}
pair<Node*, Node*> split(Node* n, int k) { // splits so left side has k nodes
	if (!n) return {};
	if (cnt(n->l) >= k) { // "n->val >= k" for lower_bound(k)
		auto pa = split(n->l, k);
		n->l = pa.second;
		n->recalc();
		return {pa.first, n};
	} else {
		auto pa = split(n->r, k - cnt(n->l) - 1); // and just "k"
		n->r = pa.first;
		n->recalc();
		return {n, pa.second};
	}
}
Node* merge(Node* l, Node* r) {
	if (!l) return r;
	if (!r) return l;
	if (l->y > r->y) {
		l->pushdown();
		l->r = merge(l->r, r);
		l->recalc();
		return l;
	} else {
		r->pushdown();
		r->l = merge(l, r->l);
		r->recalc();
		return r;
	}
}
Node* ins(Node* t, Node* n, int pos) {// insert so node is in pos (0-idx)
	auto pa = split(t, pos);
	return merge(merge(pa.first, n), pa.second);
}
Node* remove(Node* t, int pos) {
	auto pa = split(t, pos);
	auto pb = split(pa.second, 1);
	return merge(pa.first, pb.second);
}
// Example application: move the range [l, r) to index k
void move(Node*& t, int l, int r, int k) {
	Node *a, *b, *c;
	tie(a,b) = split(t, l); tie(b,c) = split(b, r - l);
	//b has the range [l, r)
	if (k <= l) t = merge(ins(a, b, k), c);
	else t = merge(a, ins(c, b, k - r));
}
