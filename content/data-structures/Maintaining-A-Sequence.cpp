#define USE_POINTER

#ifdef USE_POINTER

#include <cstdio>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cassert>
#include <ostream>
#include <vector>
#include <tuple>
#include <array>

// Sequence BBST (Implementation)
struct AVL
{
public:
	struct Node
	{
	public:
		int s, h;
		bool rev;
		std::array<Node*, 2> c;
		Node(): s(1), h(1), rev(0), c{} {}

		/* Customizable portion
		 */
		__int16_t val;
		int sum;
		std::array<int, 3> best;
		bool set;
		__int16_t set_val;
		Node(int v): Node()
		{
			val = v;
			sum = v;
			best = {v, v, v};
			set = 0;
		}

		void do_rev()
		{
			rev ^= 1;
			std::swap(c[0], c[1]);
			std::swap(best[1], best[2]);
		}
		void up()
		{
			sum = val;
			best = {val, val, val};

			if(c[0])
			{
				best = {
					std::max({best[0], c[0]->best[0], c[0]->best[2] + best[1]}),
					std::max(c[0]->best[1], c[0]->sum + best[1]),
					std::max(best[2], c[0]->best[2] + sum)
				};
				sum += c[0]->sum;
			}
			if(c[1])
			{
				best = {
					std::max({best[0], c[1]->best[0], best[2] + c[1]->best[1]}),
					std::max(best[1], sum + c[1]->best[1]),
					std::max(c[1]->best[2], c[1]->sum + best[2])
				};
				sum += c[1]->sum;
			}
		}
		void down()
		{
			if(set)
			{
				if(c[0]) c[0]->sval(set_val);
				if(c[1]) c[1]->sval(set_val);
				set = 0;
			}
		}
		void sval(__int16_t v)
		{
			val = v;
			sum = s * v;
			set = 1;
			set_val = v;
			if(v>0) best = {sum, sum, sum};
			else best = {v, v, v};
		}
		friend std::ostream& operator << (std::ostream &os, Node const& n)
		{
			os << "Node(" << n.s << ", " << n.h
				<< ", rev=" << n.rev
				<< ", val=" << n.val
				<< ", sum=" << n.sum
				<< ", set=" << n.set
				<< ", best={" << n.best[0] << ", " << n.best[1] << ", " << n.best[2] << "}"
				<< ", c={";
			if(n.c[0]) os << *n.c[0];
			os << ", ";
			if(n.c[1]) os << *n.c[1];
			return os<< "})";
		}
		/* Customizable portion (end)
		 */
	};
private:
	static int gh(Node *n) {return n ? n->h : 0;}
	static int gs(Node *n) {return n ? n->s : 0;}
	static void up(Node *n)
	{
		n->h = std::max(gh(n->c[0]), gh(n->c[1])) + 1;
		n->s = gs(n->c[0]) + gs(n->c[1]) + 1;
		// Other utilities
		n->up();
	}
	static void down(Node *n)
	{
		if(n->rev)
		{
			if(n->c[0]) n->c[0]->do_rev();
			if(n->c[1]) n->c[1]->do_rev();
			n->rev = 0;
		}
		// Other utilities
		n->down();
	}
	static Node* rotate(Node *n, int d) // Remove initial down(n) for runtime sake
	{
		Node *o = n->c[d];
		down(o); n->c[d] = o->c[!d], o->c[!d] = n;
		up(n), up(o);
		return o;
	}
	static Node* balance(Node *n) // Remove initial down(n) for runtime sake
	{
		up(n);
		int diff = gh(n->c[0]) - gh(n->c[1]);
		bool d;
		if(diff >= 2) d = 0;
		else if(diff <= -2) d = 1;
		else return n;
		down(n->c[d]);
		if(gh(n->c[d]->c[d]) + 1 < gh(n->c[d])) n->c[d] = rotate(n->c[d], !d);
		return rotate(n, d);
	}
	static Node *merge_root(Node *l, Node *n, Node *r)
	{
		if(gh(l) + 1 < gh(r))
			return down(r), r->c[0] = merge_root(l, n, r->c[0]), balance(r);
		else if(gh(r) + 1 < gh(l))
			return down(l), l->c[1] = merge_root(l->c[1], n, r), balance(l);
		else
			return n->c = {l, r}, balance(n);
	}
	static std::tuple<Node*, Node*> split(Node *n, int k)
	{
		if(n) down(n);
		if(k == 0) return {nullptr, n};
		if(k == n->s) return {n, nullptr};
		if(k <= gs(n->c[0]))
		{
			auto [l, r] = split(n->c[0], k);
			return {l, merge_root(r, n, n->c[1])};
		}
		else
		{
			auto [l, r] = split(n->c[1], k - gs(n->c[0]) - 1);
			return {merge_root(n->c[0], n, l), r};
		}
	}
	static Node* merge(Node* l, Node* r)
	{
		if(!r) return l;
		auto [x, nr] = split(r, 1);
		return merge_root(l, x, nr);
	}
public:
	Node *root;
	AVL(Node* _root=nullptr): root(_root) {}
	static AVL merge(AVL a, AVL b)
	{
		return AVL(merge(a.root, b.root));
	}
	std::tuple<AVL, AVL> split(int k)
	{
		auto [l, r] = split(root, k);
		return {AVL(l), AVL(r)};
	}
	/* Optional */
	static AVL merge_root(AVL a, AVL b, AVL c)
	{
		return AVL(merge_root(a.root, b.root, c.root));
	}
	friend std::ostream& operator << (std::ostream& os, AVL const& n)
	{
		if(n.root)
			return os << "AVL(" << *n.root << ")";
		return os << "AVL()";
	}
	/* END Optional */
};

#ifdef LOCAL
int const MV = 100;
#else
int const MV = 5e5+10;
#endif
AVL::Node pool[MV];
std::vector<int> open;

void dealloc(AVL t)
{
	auto dfs=[&](auto const& dfs, AVL::Node *n)->void
	{
		if(!n) return;
		dfs(dfs, n->c[0]);
		dfs(dfs, n->c[1]);
		open.push_back(n - pool);
	};
	dfs(dfs, t.root);
}
AVL make_avl(int v)
{
	assert(!open.empty());
	int k = open.back(); open.pop_back();
	pool[k] = AVL::Node(v);
	return pool + k;
}
AVL make_avl(std::vector<int> const& v)
{
	auto dfs=[&](auto const& dfs, int l, int r)->AVL
	{
		if(r-l==0) return AVL();
		if(r-l==1) return make_avl(v[l]);
		int m = l+(r-l)/2;
		return AVL::merge_root(dfs(dfs, l, m), make_avl(v[m]), dfs(dfs, m+1, r));
	};
	return dfs(dfs, 0, v.size());
}
int main()
{
	open.resize(MV);
	std::iota(open.begin(), open.end(), 0);
	int N, M;
	scanf("%d%d", &N, &M);
	AVL cur;
	for(int i=0;i<N;++i)
	{
		int v;
		scanf("%d", &v);
		cur = AVL::merge(cur, make_avl(v));
	}

	for(int i=0;i<M;++i)
	{
		char s[20];
		scanf(" %s", s);
		if(s[0] == 'I')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			if(t == 0) continue;
			std::vector<int> v(t);
			for(int& x: v) scanf("%d", &x);
			AVL n = make_avl(v);
			auto [l, r] = cur.split(p);
			cur = AVL::merge(l, AVL::merge(n, r));
		}
		else if(s[0] == 'D')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			dealloc(r1);
			cur = AVL::merge(l, r2);
		}
		else if(s[0] == 'M' && s[2] == 'K')
		{
			int p, t, c;
			scanf("%d%d%d", &p, &t, &c);
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			if(r1.root) r1.root->sval(c);
			cur = AVL::merge(l, AVL::merge(r1, r2));
		}
		else if(s[0] == 'R')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			if(r1.root) r1.root->do_rev();
			cur = AVL::merge(l, AVL::merge(r1, r2));
		}
		else if(s[0] == 'G')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			if(r1.root)
				printf("%d\n", r1.root->sum);
			else
				printf("0\n");
			cur = AVL::merge(l, AVL::merge(r1, r2));
		}
		else if(s[0] == 'M')
		{
			assert(cur.root);
			printf("%d\n", cur.root->best[0]);
		}
		//std::cerr << "OK " << cur << std::endl;
	}
	return 0;
}

#else

/*
 * Pointer implementation
 */

#include <cstdio>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cassert>
#include <ostream>
#include <vector>
#include <tuple>
#include <array>

template<typename T> bool ckmax(T& a, T const& b) {return b>a?a=b:1;}
template<typename T> bool ckmin(T& a, T const& b) {return b<a?a=b:1;}

struct AVL
{
public:
	struct Node
	{
	public:
		int s, h;
		std::array<int, 2> c;
		Node(): s(1), h(1), c{-1, -1} {}
		/* Customize */
		bool rev;
		int val, sum;
		bool set;
		int set_val;
		int bl, b, br;

		void do_rev()
		{
			rev ^= 1;
			std::swap(c[0], c[1]);
			std::swap(bl, br);
		}
		void do_set(int v)
		{
			val = v;
			sum = v * s;
			set = 1;
			set_val = v;
			if(v > 0) bl=sum, b=sum, br=sum;
			else bl=v, b=v, br=v;
		}
		Node(int v): Node()
		{
			rev = 0;
			val = v;
			sum = v;
			set = 0;
			bl = b = br = val;
		}
		void up()
		{
			sum = val;
			bl = b = br = val;
			if(c[0] != -1)
			{
				Node &n = N[c[0]];
				ckmax(b, std::max(n.b, n.br + bl));
				bl = std::max(n.sum + bl, n.bl);
				ckmax(br, sum + n.br);
				sum += n.sum;
			}
			if(c[1] != -1)
			{
				Node &n = N[c[1]];
				ckmax(b, std::max(n.b, br + n.bl));
				ckmax(bl, sum + n.bl);
				br = std::max(br + n.sum, n.br);
				sum += n.sum;
			}
		}
		void down()
		{
			if(rev)
			{
				if(c[0] != -1) N[c[0]].do_rev();
				if(c[1] != -1) N[c[1]].do_rev();
				rev = 0;
			}
			if(set)
			{
				if(c[0] != -1) N[c[0]].do_set(set_val);
				if(c[1] != -1) N[c[1]].do_set(set_val);
				set = 0;
			}
		}

		Node nl() const {return N[c[0]];}
		Node nr() const {return N[c[1]];}
		friend std::ostream& operator << (std::ostream &os, Node const& n)
		{
			os << "Node(" << n.s << ", " << n.h
				<< ", rev=" << n.rev
				<< ", val=" << n.val
				<< ", sum=" << n.sum
				<< ", set=" << n.set
				<< ", set_val=" << n.set_val
				<< ", best={" << n.bl << ", " << n.b << ", " << n.br << "}"
				<< ", c={";
			if(n.c[0] != -1) os << n.nl();
			os << ", ";
			if(n.c[1] != -1) os << n.nr();
			return os<< "})";
		}
		/* Customize */
	};
private:
	static std::vector<Node> N;
	static std::vector<int> bank;
	static int gh(int n) {return n != -1 ? N[n].h : 0;}
	static int gs(int n) {return n != -1 ? N[n].s : 0;}
	static void up(int n)
	{
		N[n].h = std::max(gh(N[n].c[0]), gh(N[n].c[1])) + 1;
		N[n].s = gs(N[n].c[0]) + gs(N[n].c[1]) + 1;
		N[n].up(); // Other utilities
	}
	static void down(int n)
	{
		N[n].down(); // Other utilities
	}
	static int rotate(int n, int d)
	{
		int o = N[n].c[d];
		down(o); N[n].c[d] = N[o].c[!d], N[o].c[!d] = n;
		up(n), up(o);
		return o;
	}
	static int balance(int n)
	{
		up(n);
		int diff = gh(N[n].c[0]) - gh(N[n].c[1]);
		bool d;
		if(diff >= 2) d = 0;
		else if(diff <= -2) d = 1;
		else return n;
		down(N[n].c[d]);
		if(gh(N[N[n].c[d]].c[d]) + 1 < gh(N[n].c[d])) N[n].c[d] = rotate(N[n].c[d], !d);
		return rotate(n, d);
	}
	static int merge_root(int l, int n, int r)
	{
		if(gh(l) + 1 < gh(r))
			return down(r), N[r].c[0] = merge_root(l, n, N[r].c[0]), balance(r);
		else if(gh(r) + 1 < gh(l))
			return down(l), N[l].c[1] = merge_root(N[l].c[1], n, r), balance(l);
		else
			return N[n].c = {l, r}, balance(n);
	}
	static std::tuple<int, int> split(int n, int k)
	{
		if(n != -1) down(n);
		if(k == 0) return {-1, n};
		if(k == N[n].s) return {n, -1};
		if(k <= gs(N[n].c[0]))
		{
			auto [l, r] = split(N[n].c[0], k);
			return {l, merge_root(r, n, N[n].c[1])};
		}
		else
		{
			auto [l, r] = split(N[n].c[1], k - gs(N[n].c[0]) - 1);
			return {merge_root(N[n].c[0], n, l), r};
		}
	}
	static int merge(int l, int r)
	{
		if(r == -1) return l;
		auto [x, nr] = split(r, 1);
		return merge_root(l, x, nr);
	}
	static void clear(int n)
	{
		if(n == -1)
			return;
		clear(N[n].c[0]);
		clear(N[n].c[1]);
		bank.push_back(n);
	}
	AVL(int v): root(v) {}
public:
	int root;
	AVL(): root(-1) {}
	AVL& operator += (AVL const& o) {root = merge(root, o.root); return *this;}
	friend AVL operator + (AVL a, AVL const& b) {return a += b;}
	std::tuple<AVL, AVL> split(int k)
	{
		auto [l, r] = split(root, k);
		return {AVL(l), AVL(r)};
	}
	void clear()
	{
		clear(root);
		root = -1;
	}
	Node& get_root() {return N[root];}

	/* Optional creators */
	static AVL make_avl(Node n)
	{
		int root;
		if(bank.empty()) root = N.size(), N.push_back(n);
		else root = bank.back(), bank.pop_back(), N[root] = n;
		return AVL(root);
	}
	static AVL make_avl(std::vector<Node> const& n)
	{
		auto dfs=[&](auto const& dfs, int l, int r)->int
		{
			if(r-l==0) return -1;
			if(r-l==1) return make_avl(n[l]).root;
			int m=l+(r-l)/2;
			return merge_root(dfs(dfs, l, m), make_avl(n[m]).root, dfs(dfs, m+1, r));
		};
		return AVL(dfs(dfs, 0, n.size()));
	}
	/* End optional creators */
};
typedef AVL::Node Node;
std::vector<Node> AVL::N;
std::vector<int> AVL::bank;

AVL make_avl(int v)
{
	Node n = Node(v);
	return AVL::make_avl(n);
}
AVL make_avl(std::vector<int> const& v)
{
	std::vector<Node> n(v.size());
	for(int i=0;i<v.size();++i)
		n[i] = Node(v[i]);
	return AVL::make_avl(n);
}

int main()
{
	int N, M;
	scanf("%d%d", &N, &M);
	std::vector<int> v(N);
	for(int &x: v) scanf("%d", &x);
	AVL cur = make_avl(v);
	//std::cout << cur.get_root() << '\n';
	for(int i=0;i<M;++i)
	{
		char s[100];
		scanf(" %s", s);
		if(s[0] == 'I')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			auto [l, r] = cur.split(p);
			std::vector<int> v(t);
			for(int &x: v) scanf("%d", &x);
			cur = l + make_avl(v) + r;
		}
		else if(s[0] == 'D')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			r1.clear();
			cur = l + r2;
		}
		else if(s[0] == 'M' && s[2] == 'K')
		{
			int p, t, c;
			scanf("%d%d%d", &p, &t, &c);
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			r1.get_root().do_set(c);
			cur = l + (r1 + r2);
		}
		else if(s[0] == 'R')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			r1.get_root().do_rev();
			cur = l + (r1 + r2);
		}
		else if(s[0] == 'G')
		{
			int p, t;
			scanf("%d%d", &p, &t);
			if(t == 0)
			{
				printf("0\n");
				continue;
			}
			--p;
			auto [l, r] = cur.split(p);
			auto [r1, r2] = r.split(t);
			printf("%d\n", r1.get_root().sum);
			cur = l + (r1 + r2);
		}
		else
			printf("%d\n", cur.get_root().b);
		//std::cout << cur.get_root() << '\n';
	}
	return 0;
}
#endif
