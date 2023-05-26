/**
 * Author: Jacob
 * Date: back in the day
 * License: CC0
 * Source: me
 * Description: 2D point-update range-query segtree supporting $10^9$ coordinates (IOI 2013 game)
 * Time: $O(\log^2 N)$
 * Status: tested on IOI
 */

#pragma once

#define DEFAULT 0ll
ll func(ll a, ll b) {return max(a, b);} // associative func
struct SegTree2D {
	int R,C,root;
	vector<int> l,r,b,e,st;
	vector<ll> v;
	inline int mid(int x,int y) {return ((x+y)>>1);}
	SegTree2D(int _R,int _C):R(_R),C(_C) {
		l.pb(0),r.pb(0),b.pb(0),e.pb(0),st.pb(0),v.pb(DEFAULT);
		R=_R;
		root=alloc2(0,R-1);
	}
	int alloc(int _b,int _e,ll _v) {
		l.pb(0),r.pb(0),b.pb(_b),e.pb(_e),st.pb(0),v.pb(_v);
		return sz(b)-1;
	}
	void lca(int b,int e,int ob,int oe,int i,int &nb,int &ne) {
		int m=mid(b,e);
		if ((i<=m&&ob>m)||(i>m&&oe<=m)) nb=b,ne=e;
		else (i>m)?lca(m+1,e,ob,oe,i,nb,ne):lca(b,m,ob,oe,i,nb,ne);
	}
	void up(int x) {v[x]=func(v[l[x]],v[r[x]]);}
	void update1(int x,int i,ll nv) {
		if (b[x]>i||e[x]<i) return;
		if (b[x]==e[x]) {
			v[x]=nv;
			return;
		}
		int m=mid(b[x],e[x]);
		if (i<=m) {
			if (l[x]) {
				if (b[l[x]]<=i && i<=e[l[x]]) update1(l[x],i,nv);
				else {
					int nb,ne;
					lca(0,C-1,b[l[x]],e[l[x]],i,nb,ne);
					int y=l[x];
					l[x]=alloc(nb,ne,DEFAULT);
					if (i>mid(nb,ne)) l[l[x]]=y,r[l[x]]=alloc(i,i,nv);
					else r[l[x]]=y,l[l[x]]=alloc(i,i,nv);
					up(l[x]);
				}
			} else l[x]=alloc(i,i,nv);
		} else {
			if (r[x]) {
				if (b[r[x]]<=i && i<=e[r[x]]) update1(r[x],i,nv);
				else {
					int nb,ne;
					lca(0,C-1,b[r[x]],e[r[x]],i,nb,ne);
					int y=r[x];
					r[x]=alloc(nb,ne,DEFAULT);
					if (i>mid(nb,ne)) l[r[x]]=y,r[r[x]]=alloc(i,i,nv);
					else r[r[x]]=y,l[r[x]]=alloc(i,i,nv);
					up(r[x]);
				}
			} else r[x]=alloc(i,i,nv);
		}
		up(x);
	}
	ll query1(int x,int qb,int qe) {
		if (!x) return DEFAULT;
		if (b[x]>qe||e[x]<qb) return DEFAULT;
		if (b[x]>=qb&&e[x]<=qe) return v[x];
		return func(query1(l[x],qb,qe),query1(r[x],qb,qe));
	}
	int alloc2(int _b,int _e) {
		int newnode = alloc(0,C-1,DEFAULT);
		l.pb(0),r.pb(0),b.pb(_b),e.pb(_e),v.pb(DEFAULT),st.pb(newnode);
		return sz(b)-1;
	}
	void update2(int x,int i,int j,ll nv) {
		if (b[x]>i||e[x]<i) return;
		if (b[x]==e[x]) update1(st[x],j,nv);
		else {
			int m=mid(b[x],e[x]);
			if (!l[x]) {
				l[x]=alloc2(b[x],m);
				r[x]=alloc2(m+1,e[x]);
			}
			if (i<=m) update2(l[x],i,j,nv);
			else update2(r[x],i,j,nv);
			update1(st[x],j,func(query1(st[l[x]],j,j),query1(st[r[x]],j,j)));
		}
	}
	ll query2(int x,int rb,int re,int cb,int ce) {
		if (!x) return DEFAULT;
		if (b[x]>re||e[x]<rb) return DEFAULT;
		if (b[x]>=rb&&e[x]<=re) return query1(st[x],cb,ce);
		return func(query2(l[x],rb,re,cb,ce),query2(r[x],rb,re,cb,ce));
	}
	void update(int p,int q,ll k) {update2(root,p,q,k);}
	ll query(int p,int q,int u,int v) {return query2(root,p,u,q,v);}
};