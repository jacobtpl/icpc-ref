/**
 * Author: Benq
 * Description: Full solution to "Dynamic Tree Test": link-cut tree that also
 *  keeps virtual subtrees (c[2], c[3] and c[4]), so it supports path and subtree
 *  set/add and min/max/sum queries, reroot and change-parent.
 *  INFO[0] is the path aggregate, INFO[1] the aggregate of everything hanging
 *  off the path. A whole program, not a header: adapt main().
 * Time: O(log N) amortized per operation, large constant
 * Status: checked against a brute force on small random cases (test not in repo)
 */
#include <bits/stdc++.h>
using namespace std;
typedef vector<int> vi;
typedef pair<int,int> pi;

#define FOR(i, a, b) for (int i = (a); i < (b); i++)
#define F0R(i, a) for (int i = 0; i < (a); i++)
#define FORd(i,a,b) for (int i = (b)-1; i >= (a); i--)
#define F0Rd(i,a) for (int i = (a)-1; i >= 0; i--)
#define each(a, x) for (auto& a : x)

#define mp make_pair
#define pb push_back
#define f first
#define s second
#define lb lower_bound
#define ub upper_bound

#define sz(x) (int)x.size()
#define beg(x) x.begin()
#define en(x) x.end()
#define all(x) beg(x), en(x)
#define resz resize

template<class T> void ckmin(T &a, T b) { a = min(a, b); }
template<class T> void ckmax(T &a, T b) { a = max(a, b); }

typedef struct snode* sn;

struct info {
    int sz, mn, mx; long long sum;
    info (int v) {
        if (v == INT_MAX) {
            sz = 0; sum = 0;
            mn = INT_MAX, mx = INT_MIN;
        } else {
            sz = 1; sum = mn = mx = v;
        }
    }
    info() : info(INT_MAX) {}
};

info& operator+=(info& a, const info& b) {
    a.sz += b.sz, a.sum += b.sum;
    ckmin(a.mn,b.mn); ckmax(a.mx,b.mx);
    return a;
}

struct snode {
    //////// VARIABLES
    int id, val; // value in node
    sn p, c[5]; // parent, children
    bool flip = 0; 
    info INFO[2]; 
    int nextNum[2], lazy[2];
    
    snode(int ID, int v) {
        id = ID; val = v; 
        p = NULL; F0R(i,5) c[i] = NULL;
        nextNum[0] = nextNum[1] = INT_MAX; lazy[0] = lazy[1] = 0;
        calc();
    }
    
    //////// SPLAY TREE OPERATIONS
    void prop() {
        if (flip) {
            swap(c[0],c[1]); 
            F0R(i,2) if (c[i]) c[i]->flip ^= 1;
            flip = 0; 
        }
        if (nextNum[1] != INT_MAX) { // modify everything in subtree EXCEPT 01 paths 
            if (INFO[1].sz) {
                INFO[1].sum = (long long)nextNum[1]*INFO[1].sz; 
                INFO[1].mn = INFO[1].mx = nextNum[1];
            }
            F0R(i,5) if (c[i]) {
                c[i]->nextNum[1] = nextNum[1], c[i]->lazy[1] = 0;
                if (i >= 2) c[i]->nextNum[0] = nextNum[1], c[i]->lazy[0] = 0;
            }
            nextNum[1] = INT_MAX;
        }
        if (lazy[1] != 0) {
            if (INFO[1].sz) {
                INFO[1].sum += (long long)lazy[1]*INFO[1].sz;
                INFO[1].mn += lazy[1], INFO[1].mx += lazy[1];
            }
            F0R(i,5) if (c[i]) {
                c[i]->lazy[1] += lazy[1];
                if (i >= 2) c[i]->lazy[0] += lazy[1];
            }
            lazy[1] = 0;
        }
        if (nextNum[0] != INT_MAX) { // modify ONLY 01 paths 
            val = nextNum[0];
            INFO[0].sum = (long long)nextNum[0]*INFO[0].sz;
            INFO[0].mn = INFO[0].mx = nextNum[0];
            F0R(i,2) if (c[i]) c[i]->nextNum[0] = nextNum[0], c[i]->lazy[0] = 0;
            nextNum[0] = INT_MAX;
        }
        if (lazy[0] != 0) {
            val += lazy[0];
            INFO[0].sum += (long long)lazy[0]*INFO[0].sz;
            INFO[0].mn += lazy[0], INFO[0].mx += lazy[0];
            F0R(i,2) if (c[i]) c[i]->lazy[0] += lazy[0];
            lazy[0] = 0;
        }
    }
    void calc() { 
        F0R(i,5) if (c[i]) c[i]->prop();
        INFO[0] = info(val); INFO[1] = info(INT_MAX); 
        F0R(i,5) if (c[i]) {
            INFO[1] += c[i]->INFO[1];
            if (i >= 2) INFO[1] += c[i]->INFO[0];
            else INFO[0] += c[i]->INFO[0];
        }
    }
    int dir() {
        if (!p) return -2;
        F0R(i,5) if (p->c[i] == this) return i;
        assert(0);
    }
    bool isRoot() { int d = dir(); return d == -2 || d == 4; }
    
    friend void setLink(sn x, sn y, int d) {
        if (y) y->p = x;
        if (d >= 0) x->c[d] = y;
    }
    void rot() { // assume p and p->p propagated
        assert(!isRoot()); int x = dir(); sn pa = p;
        setLink(pa->p, this, pa->dir());
        setLink(pa, c[x^1], x);
        setLink(this, pa, x^1);
        pa->calc(); calc();
    }
    bool ok0() {
        int d = dir(); return d == 0 || d == 1;
    }
    void splay() {
        while (ok0() && p->ok0() && p->p->ok0()) {
            p->p->prop(), p->prop(), prop();
            dir() == p->dir() ? p->rot() : rot();
            rot();
        }
        if (ok0() && p->ok0()) p->prop(), prop(), rot();
        if (ok0()) {
            p->prop(), prop();
            auto a = p->p, b = p->c[2], c = p->c[3]; int d = p->dir();
            p->p = p->c[2] = p->c[3] = NULL; p->calc(); rot(); 
            setLink(this,b,2); setLink(this,c,3); setLink(a,this,d); calc();
        }
        while (!isRoot() && !p->isRoot()) {
            p->p->prop(), p->prop(), prop();
            dir() == p->dir() ? p->rot() : rot();
            rot();
        }
        if (!isRoot()) p->prop(), prop(), rot();
        prop();
    }

    sn splayRight() {
        prop();
        if (!c[3]) { splay(); return this; }
        return c[3]->splayRight();
    }
    
    friend sn join(sn a, sn b) {
        if (!a) return b;
        a->splay(); a = a->splayRight();
        setLink(a,b,3); a->calc(); return a;
    }
    
    //////// LINK CUT TREE BASE OPERATIONS
    void access() { // bring this to top of tree, left subtree of this is now path to root
        int it = 0;
        for (sn v = this, pre = NULL; v; v = v->p) {
            it ++;
            v->splay(); auto c = v->c[1]; 
            if (c) assert(!c->c[2] && !c->c[3]);
            if (pre) pre->prop();
            if (pre) {
                assert(v->c[4] == pre);
                auto a = pre->c[2], b = pre->c[3]; 
                if (a) a->p = NULL;
                if (b) b->p = NULL;
                pre->c[2] = pre->c[3] = NULL; pre->calc();
                if (c) c->p = NULL;
                setLink(v,join(join(a,b),c),4); 
            } else {
                if (c) c->p = NULL;
                if (v->c[4]) v->c[4]->p = NULL;
                setLink(v,join(c,v->c[4]),4);
            }
            v->c[1] = pre; v->calc(); pre = v;
        }
        splay(); assert(!c[1]);
    }
    void makeRoot() { access(); flip ^= 1; }
    
    //////// LINK CUT TREE QUERIES
    friend sn lca(sn x, sn y) {
        if (x == y) return x;
        x->access(); y->access(); 
        if (!x->p) return NULL;
        x->splay(); return x->p ? x->p : x;
    }
    friend bool connected(sn x, sn y) { return lca(x,y); }
    
    //////// LINK CUT TREE MODIFICATIONS
    friend bool link(sn x, sn y) { // make y parent of x
        if (connected(x,y)) exit(2);
        x->makeRoot(); setLink(y,join(x,y->c[4]),4); y->calc();
        return 1;
    }
    friend bool cut(sn x, sn y) {
        x->makeRoot(); y->access(); 
        if (y->c[0] != x || x->c[0] || x->c[1]) exit(3);
        x->p = y->c[0] = NULL; y->calc(); 
        return 1;// ??
    }
    
    void propAll() {
        prop();
        F0R(i,5) if (c[i]) c[i]->propAll();
    }
};

sn getPar(sn x) {
    x->access(); x = x->c[0];
    while (1) {
        x->prop();
        if (!x->c[1]) return x;
        x = x->c[1];
    }
    return x;
}

int N,M,root;
const int MX = 100005;
sn LCT[MX];

void pathInc(int y, int z) {
    LCT[y]->access(); LCT[y]->lazy[0] = z;
}

void pathMod(int y, int z) {
    LCT[y]->access(); LCT[y]->nextNum[0] = z;
}

void init() {
    cin.tie(0)->sync_with_stdio(0);
    cin >> N >> M;
    vector<pi> ed(N-1);
    for (auto& e : ed) cin >> e.f >> e.s;
    FOR(i,1,N+1) {
        int x; cin >> x;
        LCT[i] = new snode(i,x);
    }
    F0R(i,N-1) {
        int x = ed[i].f, y = ed[i].s;
        link(LCT[x],LCT[y]);
    }
    cin >> root; 
}

void PS(long long x) { cout << x << "\n"; }

int main() {
    init();
    F0R(i,M) {
        int K; cin >> K;
        if (K == 0 || K == 5) {
            int x,y; cin >> x >> y; LCT[root]->makeRoot(); LCT[x]->access();
            auto c = LCT[x]->c[4];
            if (K == 0) {
                if (c) c->prop(), c->nextNum[0] = c->nextNum[1] = y;
                LCT[x]->val = y; 
                
            } else {
                if (c) c->prop(), c->lazy[0] = c->lazy[1] = y;
                LCT[x]->val += y;
            }
            LCT[x]->calc();
        } else if (K == 1) {
            int x; cin >> x; root = x;
        } else if (K == 2 || K == 6) {
            int x,y,z; cin >> x >> y >> z; LCT[x]->makeRoot(); 
            if (K == 2) pathMod(y,z);
            else pathInc(y,z);
        } else if (K == 3 || K == 4 || K == 11) {
            int x; cin >> x; LCT[root]->makeRoot(); LCT[x]->access();
            auto ans = info(LCT[x]->val); sn c = LCT[x]->c[4];
            if (c) c->prop(), ans += c->INFO[0], ans += c->INFO[1];
            if (K == 3) PS(ans.mn);
            else if (K == 4) PS(ans.mx);
            else if (K == 11) PS(ans.sum);
        } else if (K == 7 || K == 8 || K == 10) {
            int x,y; cin >> x >> y; LCT[x]->makeRoot(); LCT[y]->access();
            auto ans = LCT[y]->INFO[0];
            if (K == 7) PS(ans.mn);
            else if (K == 8) PS(ans.mx);
            else PS(ans.sum);
        } else if (K == 9) {
            int x,y; cin >> x >> y;
            LCT[root]->makeRoot();
            if (lca(LCT[x],LCT[y]) == LCT[x]) {
                continue;
            }
            cut(getPar(LCT[x]),LCT[x]);
            link(LCT[x],LCT[y]);
        } else exit(5);
    }
}
