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
		void do_rev()
		{
			rev ^= 1;
			std::swap(c[0], c[1]);
		}
		void up() {}
		void down() {}
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
		n->up(); // Other utilities
	}
	static void down(Node *n)
	{
		if(n->rev)
		{
			if(n->c[0]) n->c[0]->do_rev();
			if(n->c[1]) n->c[1]->do_rev();
			n->rev = 0;
		}
		n->down(); // Other utilities
	}
	static Node* rotate(Node *n, int d)
	{
		Node *o = n->c[d];
		down(o); n->c[d] = o->c[!d], o->c[!d] = n;
		up(n), up(o);
		return o;
	}
	static Node* balance(Node *n)
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
};
