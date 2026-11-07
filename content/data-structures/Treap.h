/**
 * Author: Benq, Jacob Teo
 * Source: folklore
 * Description: Treap, a short self-balancing tree. Acts as an
 *  ordered set (split/ins/del by key) or as a sequential
 *  container (splitsz/inspos/delpos by 0-indexed position)
 *  with log-time split/merge. Augment by adding fields and
 *  updating calc (pull) and prop (lazy push); sz, sum and
 *  flip are examples. Nodes are never freed.
 * Time: $O(\log N)$
 * Status: stress-tested
 */
#pragma once

using pt = struct Node*;
struct Node {
	int pri, val; pt c[2]; // essential
	int sz = 1; ll sum; // for range queries
	bool flip = 0; // lazy update
	Node(int v) : pri(rand()), val(v), c{}, sum(v) {}
};
int getsz(pt x) { return x ? x->sz : 0; }
ll getsum(pt x) { return x ? x->sum : 0; }
pt prop(pt x) { // lazy propagation
	if (!x || !x->flip) return x;
	swap(x->c[0], x->c[1]);
	x->flip = 0; rep(i,0,2) if (x->c[i]) x->c[i]->flip ^= 1;
	return x;
}
pt calc(pt x) { // pull from (propagated) children
	pt a = x->c[0], b = x->c[1];
	assert(!x->flip); prop(a), prop(b);
	x->sz = 1 + getsz(a) + getsz(b);
	x->sum = x->val + getsum(a) + getsum(b);
	return x;
}
template<class F> void each(pt x, F f) { // inorder
	if (!x) return;
	prop(x); each(x->c[0], f); f(x->val); each(x->c[1], f);
}
pair<pt,pt> split(pt t, int v) { // >= v goes to the right
	if (!t) return {t,t};
	prop(t);
	if (t->val >= v) {
		auto p = split(t->c[0], v); t->c[0] = p.second;
		return {p.first, calc(t)};
	} else {
		auto p = split(t->c[1], v); t->c[1] = p.first;
		return {calc(t), p.second};
	}
}
pair<pt,pt> splitsz(pt t, int sz) { // sz nodes go to the left
	if (!t) return {t,t};
	prop(t);
	if (getsz(t->c[0]) >= sz) {
		auto p = splitsz(t->c[0], sz); t->c[0] = p.second;
		return {p.first, calc(t)};
	} else {
		auto p = splitsz(t->c[1], sz-getsz(t->c[0])-1);
		t->c[1] = p.first;
		return {calc(t), p.second};
	}
}
pt merge(pt l, pt r) { // keys in l < keys in r
	if (!l || !r) return l ? l : r;
	prop(l), prop(r); pt t;
	if (l->pri > r->pri) l->c[1] = merge(l->c[1], r), t = l;
	else r->c[0] = merge(l, r->c[0]), t = r;
	return calc(t);
}
pt ins(pt x, int v) { // insert v (set: existing v dropped)
	auto a = split(x, v), b = split(a.second, v+1);
	return merge(a.first, merge(new Node(v), b.second)); }
pt del(pt x, int v) { // delete v
	auto a = split(x, v), b = split(a.second, v+1);
	return merge(a.first, b.second); }
pt inspos(pt t, pt n, int pos) { // n becomes t[pos]
	auto a = splitsz(t, pos);
	return merge(merge(a.first, n), a.second); }
pt delpos(pt t, int pos) {
	auto a = splitsz(t, pos), b = splitsz(a.second, 1);
	return merge(a.first, b.second); }
pt rev(pt t, int l, int r) { // reverse [l, r)
	auto a = splitsz(t, l), b = splitsz(a.second, r-l);
	if (b.first) b.first->flip ^= 1;
	return merge(a.first, merge(b.first, b.second)); }
