/**
 * Author: Benq
 * Source: folklore
 * Description: A short self-balancing tree. It acts as a
 *  sequential container with log-time splits/joins, and
 *  is easy to augment with additional data. 0-index.
 * Time: $O(\log N)$
 * Status: stress-tested
 */
using pt = struct Node*;
struct Node {
	int pri, val; pt c[2]; // essential
	int sz; ll sum; // for range queries
	bool flip = 0; // lazy update
	Node(int _val) {
		pri = rand(); sum = val = _val; 
		sz = 1; c[0] = c[1] = nullptr;
	}
	~Node() { rep(i,0,2) delete c[i]; } /// assume no sharing of data
};
int getsz(pt x) { return x?x->sz:0; }
ll getsum(pt x) { return x?x->sum:0; }
pt prop(pt x) { // lazy propagation
	if (!x || !x->flip) return x;	
	swap(x->c[0],x->c[1]);
	x->flip = 0; rep(i,0,2) if (x->c[i]) x->c[i]->flip ^= 1;
	return x;
}
pt calc(pt x) {
	pt a = x->c[0], b = x->c[1];
	assert(!x->flip); prop(a), prop(b);
	x->sz = 1+getsz(a)+getsz(b);
	x->sum = x->val+getsum(a)+getsum(b);
	return x;
}
void tour(pt x, vi& v) { // print values of nodes, 
	if (!x) return; // inorder traversal
	prop(x); tour(x->c[0],v); v.pb(x->val); tour(x->c[1],v);
}
pair<pt,pt> split(pt t, int v) { // >= v goes to the right
	if (!t) return {t,t};
	prop(t);
	if (t->val >= v) {
		auto p = split(t->c[0], v); t->c[0] = p.second;
		return {p.first,calc(t)};
	} else {
		auto p = split(t->c[1], v); t->c[1] = p.first;
		return {calc(t),p.second};
	}
}
pair<pt,pt> splitsz(pt t, int sz) { // sz nodes go to left
	if (!t) return {t,t};
	prop(t);
	if (getsz(t->c[0]) >= sz) {
		auto p = splitsz(t->c[0],sz); t->c[0] = p.second;
		return {p.first,calc(t)};
	} else {
		auto p=splitsz(t->c[1],sz-getsz(t->c[0])-1); t->c[1]=p.first;
		return {calc(t),p.second};
	}
}
pt merge(pt l, pt r) { //  keys in l < keys in r
	if (!l || !r) return l?:r;
	prop(l), prop(r); pt t;
	if (l->pri > r->pri) l->c[1] = merge(l->c[1],r), t = l;
	else r->c[0] = merge(l,r->c[0]), t = r;
	return calc(t);
}
pt ins(pt x, int v) { // insert v
	auto a = split(x,v), b = split(a.second,v+1);
	return merge(a.first,merge(new Node(v),b.second)); }
pt del(pt x, int v) { // delete v
	auto a = split(x,v), b = split(a.second,v+1);
	return merge(a.first,b.second); }
pt inspos(pt t, pt n, int pos) {// insert so node is in pos (0-indexed)
	auto pa = splitsz(t, pos);
	return merge(merge(pa.first, n), pa.second); }
pt delpos(pt t, int pos) {
	auto pa = splitsz(t, pos);
	auto pb = splitsz(pa.second, 1);
	return merge(pa.first, pb.second); }