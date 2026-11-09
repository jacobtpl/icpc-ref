/**
 * Author: Benq, Jacob Teo
 * Source: folklore
 * Description: Treap, a short self-balancing tree. Acts as an
 *  ordered set (split/ins/del by key) or as a sequential
 *  container (splitsz/inspos/delpos by 0-indexed position)
 *  with log-time split/merge. Nodes are never freed.
 *  To add an aggregate (like sz, sum, mn): add a field and a
 *  line in calc. To add a lazy (like add, flip): add a tag
 *  field and an apply function like upd/rev that updates x's
 *  own fields and tag, then call it on the children in prop.
 *  Only touch tags and aggregates through those functions.
 * Usage: pt t = build({1, 2, 3}); // O(n)
 *  t = onRange(t, l, r, [&](pt x) { rev(x); }); // reverse
 *  t = onRange(t, l, r, [&](pt x) { upd(x, 5); }); // add 5
 *  t = onRange(t, l, r, [&](pt x) { ans = x->mn; }); // query
 * Time: $O(\log N)$
 * Status: stress-tested
 */
#pragma once

using pt = struct Node*;
pt calc(pt x);
struct Node {
	int pri; ll val; pt c[2]; // essential
	int sz; ll sum, mn; // aggregates, set by calc
	bool flip = 0; ll add = 0; // lazy tags, pending for children
	Node(ll v) : pri(rand()), val(v), c{} { calc(this); }
};
int getsz(pt x) { return x ? x->sz : 0; }
ll getsum(pt x) { return x ? x->sum : 0; }
ll getmn(pt x) { return x ? x->mn : LLONG_MAX; }
void upd(pt x, ll v) { // add v to every value in x's subtree
	x->val += v; x->mn += v; x->sum += v * x->sz; x->add += v;
}
void rev(pt x) { swap(x->c[0], x->c[1]); x->flip ^= 1; }
pt prop(pt x) { // push tags to children (order can matter)
	if (!x) return x;
	rep(i,0,2) if (pt y = x->c[i]) {
		if (x->flip) rev(y);
		if (x->add) upd(y, x->add);
	}
	x->flip = 0, x->add = 0;
	return x;
}
pt calc(pt x) { // pull aggregates from children
	pt a = x->c[0], b = x->c[1];
	assert(!x->flip && !x->add); prop(a), prop(b);
	x->sz = 1 + getsz(a) + getsz(b);
	x->sum = x->val + getsum(a) + getsum(b);
	x->mn = min({x->val, getmn(a), getmn(b)});
	return x;
}
template<class F> void each(pt x, F f) { // inorder, f(pt)
	if (!x) return;
	prop(x); each(x->c[0], f); f(x); each(x->c[1], f);
}
pair<pt,pt> split(pt t, ll v) { // >= v goes to the right
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
pt ins(pt x, ll v) { // set insert (an existing v is dropped)
	auto a = split(x, v), b = split(a.second, v+1);
	return merge(a.first, merge(new Node(v), b.second)); }
pt del(pt x, ll v) { // delete v
	auto a = split(x, v), b = split(a.second, v+1);
	return merge(a.first, b.second); }
pt inspos(pt t, pt n, int pos) { // n becomes t[pos]
	auto a = splitsz(t, pos);
	return merge(merge(a.first, n), a.second); }
pt delpos(pt t, int pos) {
	auto a = splitsz(t, pos), b = splitsz(a.second, 1);
	return merge(a.first, b.second); }
template<class F> pt onRange(pt t, int l, int r, F f) {
	auto a = splitsz(t, l), b = splitsz(a.second, r-l);
	if (b.first) f(b.first); // f sees the subtree for [l, r)
	return merge(a.first, merge(b.first, b.second)); }
pt build(const vector<ll>& v) { // O(n), in order
	vector<pt> st;
	for (ll x : v) {
		pt n = new Node(x), last = 0;
		while (!st.empty() && st.back()->pri < n->pri)
			last = calc(st.back()), st.pop_back();
		n->c[0] = last;
		if (!st.empty()) st.back()->c[1] = n;
		st.push_back(n);
	}
	pt root = st.empty() ? 0 : st[0];
	while (!st.empty()) calc(st.back()), st.pop_back();
	return root;
}
