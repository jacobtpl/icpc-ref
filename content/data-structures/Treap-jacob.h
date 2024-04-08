#include <cstdio>
#include <cstdlib>
#include <utility>
#include <algorithm>
using namespace std;
#define ll long long
#define mp make_pair
#define sz(x) (int)(x).size()
#define all(x) (x).begin(),(x).end()
#define pb push_back
#define ii pair<int,int>
#define INF 1000000000
#define INFLL 1000000000000000010ll
#define M 1000000007ll
#define UQ(x) (x).resize(distance((x).begin(),unique(all((x)))))
#define mid(x,y) (((x)+(y))>>1)
#define POOL 3000005
struct treap {
	int c[POOL][2],sub[POOL],rev[POOL];
	int v[POOL],lazy[POOL],mv[POOL],cur,pri[POOL];
	void init() {
		cur=1;
		sub[0]=0;
	}
	bool suck(int a, int b) {
		return rand() * (sub[a]+sub[b]) < sub[a] * RAND_MAX;
	}
	ll bigrand() {
		return rand()*RAND_MAX + rand();
	}
	void alloc(int &x,ll _v) {
		c[cur][0]=c[cur][1]=0;
		sub[cur]=1;
		rev[cur]=0;
		lazy[cur]=0ll;
		pri[cur]=bigrand();
		v[cur]=mv[cur]=_v;
		x=cur++;
	}
	void pushup(int x) {
		sub[x]=sub[c[x][0]]+sub[c[x][1]]+1;
		mv[x]=v[x];
		if (c[x][0]) mv[x]=min(mv[x],mv[c[x][0]]);
		if (c[x][1]) mv[x]=min(mv[x],mv[c[x][1]]);
	}
	void pushdown(int x) {
		if (rev[x]) {
			if (c[x][0]) 
				rev[c[x][0]]^=1,swap(c[c[x][0]][0],c[c[x][0]][1]);
			if (c[x][1]) 
				rev[c[x][1]]^=1,swap(c[c[x][1]][0],c[c[x][1]][1]);
			rev[x]=0;
		}
		if (lazy[x]) {
			if (c[x][0]) {
				lazy[c[x][0]]+=lazy[x];
				mv[c[x][0]]+=lazy[x];
				v[c[x][0]]+=lazy[x];
			}
			if (c[x][1]) {
				lazy[c[x][1]]+=lazy[x];
				mv[c[x][1]]+=lazy[x];
				v[c[x][1]]+=lazy[x];
			}
			lazy[x]=0ll;
		}
	}
	void build(int *a,int &z,int b,int e) {
		if (e<b) {
			z=0;
			return;
		}
		if (b==e) {
			alloc(z,a[b]);
			return;
		}
		int m=mid(b,e),x,y;
		build(a,x,b,m);
		build(a,y,m+1,e);
		merge(z,x,y);
		pushup(z);
	}
	void merge(int &z,int x,int y) {
		if (!x || !y) {
			z=x|y;
			return;
		}
		if (pri[x]<pri[y]) {
			pushdown(x);
			merge(c[x][1],c[x][1],y);
			pushup(z=x);
		} else {
			pushdown(y);
			merge(c[y][0],x,c[y][0]);
			pushup(z=y);
		}
	}
	void split(int z,int &x,int &y,int s) {
		if (s==0) {
			y=z;
			x=0;
			return;
		}
		pushdown(z);
		if (sub[c[z][0]]>=s) {
			y=z;
			split(c[z][0],x,c[y][0],s);
			pushup(y);
		} else {
			x=z;
			split(c[z][1],c[x][1],y,s-sub[c[z][0]]-1);
			pushup(x);
		}
	}
	void output(int x) {
		if (c[x][0]) output(c[x][0]);
		printf("(%lld, %lld to %lld %lld) ", 
				v[x],mv[x],v[c[x][0]],v[c[x][1]]);
		if (c[x][1]) output(c[x][1]);
	}
};
treap t;
int n,q,o[100005];
char s[10];
void work() {
	t.init();
	for (int i=0;i<n;i++) {
		scanf("%d",&o[i]);
	}
	int root;
	t.build(o,root,0,n-1);
	int a,b,c,x,y,z;
	scanf("%d",&q);
	for (int kk=0;kk<q;kk++) {
		scanf("%s",s);
		if (s[0]=='A') {
			scanf("%d%d%d",&a,&b,&c);
			t.split(root,x,y,a-1);
			t.split(y,y,z,b-a+1);
			t.lazy[y]+=c;
			t.v[y]+=c;
			t.mv[y]+=c;
			t.merge(x,x,y);
			t.merge(root,x,z);
		} else if (s[0]=='I') {
			scanf("%d%d",&a,&c);
			t.split(root,x,z,a);
			t.alloc(y,c);
			t.merge(x,x,y);
			t.merge(root,x,z);
		} else if (s[0]=='D') {
			scanf("%d",&a);
			t.split(root,x,y,a-1);
			t.split(y,y,z,1);
			t.merge(root,x,z);
		} else if (s[0]=='M') {
			scanf("%d%d",&a,&b);
			t.split(root,x,y,a-1);
			t.split(y,y,z,b-a+1);
			printf("%d\n", t.mv[y]);
			t.merge(x,x,y);
			t.merge(root,x,z);
		} else if (s[0]=='R' && s[3]=='E') {
			scanf("%d%d",&a,&b);
			t.split(root,x,y,a-1);
			t.split(y,y,z,b-a+1);
			t.rev[y]^=1;
			swap(t.c[y][0],t.c[y][1]);
			t.merge(x,x,y);
			t.merge(root,x,z);
		} else if (s[0]=='R' && s[3]=='O') {
			scanf("%d%d%d",&a,&b,&c);
			int l=b-a+1;
			c%=l;
			c=(c+l)%l;
			if (!c) continue;
			t.split(root,x,y,a-1);
			t.split(y,y,z,b-a+1);
			int u;
			t.split(y,y,u,l-c);
			t.merge(y,u,y);
			t.merge(y,y,z);
			t.merge(root,x,y);
		}
	}
}
int main() {
	srand(89302);
	while(scanf("%d",&n)!=EOF) {
		work();
	}
}
