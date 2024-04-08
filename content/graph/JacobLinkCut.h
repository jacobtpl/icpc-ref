/**
 * Author: Jacob Teo
 * Date: 2024
 * Description: Link-cut tree with evert, node update and path sum.
 * Time: All operations take amortized O(\log N).
 */
#pragma once

struct node {
	node *p,*c[2];
	ll v, subv;
	int rev;
	node():p(NULL),rev(0),v(0),subv(0){c[0]=c[1]=NULL;}
	node(node *_p):p(_p),rev(0),v(0),subv(0){c[0]=c[1]=NULL;}
	int state() {
		if (!p) return -1;
		if (this==p->c[0]) return 0;
		else if (this==p->c[1]) return 1;
		return -1;
	}
	bool isroot() { return state()==-1; }
	node* prop() { // propagate rev, pushdown
		if (rev) {
			rev=0;
			swap(c[0],c[1]);
			if (c[0]) c[0]->rev^=1;
			if (c[1]) c[1]->rev^=1;
		}
		return this;
	}
	void set() { // update any subtrees, pullup
		subv=v;
		if (c[0]) subv+=c[0]->subv;
		if (c[1]) subv+=c[1]->subv;
	}
	void res() {
		if (p->p) p->p->prop();
		p->prop(); prop();
	}
	friend void linknode(node *px, node *x, int d) {
		if (px && d!=-1) px->c[d] = x;
		if (x) x->p = px;
	}
	void rot() {
		int d=state(),d2=p->state();
		node *b=c[!d],*pp=p,*gp=p->p;
		linknode(pp,b,d);
		linknode(this,pp,!d);
		linknode(gp,this,d2);
		c[!d]->set(); set();
	}
	node* splay() {
		while(!isroot()) {
			res();
			if (p->isroot()) rot();
			else if (state()==p->state()) p->rot(),rot();
			else rot(),rot();
		}
		return prop();
	}
	node* find_min() {
		node *x=this;
		while(x->prop()->c[0]) x=x->c[0];
		return x->splay();
	}
};
void access(node *x) {
	node *prev=NULL;
	for (node *z=x;z;z=z->p) {
		z->splay();
		z->c[1]=prev;
		z->set();
		prev=z;
	}
	x->splay();
}
void evert(node *x) {access(x); x->rev^=1;}
void link(node *x,node *y) {
	evert(x); access(x); access(y);
	y->c[1]=x; x->p=y; y->set();
}
void cut(node *x,node *y) {
	evert(y); access(x); access(y);
	x->p=NULL;
}
node* find_root(node *x) {access(x); return x->find_min();}
void update(node *x, ll v) {access(x); x->v += v; x->splay();}
ll query(node *x, node *y) {evert(x); access(y); return y->subv;}
