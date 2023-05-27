/**
 * Author: Siyong
 * Date: 2023-05-25
 * License: CC0
 * Source: me
 * Description: AVL tree with split (maintaining a sequence)
 * Time: O(\log N)
 */
struct AVL {
public:
	struct Node {
	public:
		int s, h; /* Customize */
		std::array<int, 2> c;
		Node() : s(1), h(1), c{ -1, -1 } {}
		void up(); /* Customize */
		void down();
	};
private:
	static std::vector<Node> N;
	static std::vector<int> bank;
	static int gh(int n) { return n != -1 ? N[n].h : 0; }
	static int gs(int n) { return n != -1 ? N[n].s : 0; }
	static void up(int n) {
		N[n].h = std::max(gh(N[n].c[0]), gh(N[n].c[1])) + 1;
		N[n].s = gs(N[n].c[0]) + gs(N[n].c[1]) + 1;
		N[n].up();
	}
	static void down(int n) { N[n].down(); }
	static int rotate(int n, int d) {
		int o = N[n].c[d]; down(o);
		N[n].c[d] = N[o].c[!d], N[o].c[!d] = n;
		up(n), up(o);
		return o;
	}
	static int balance(int n) {
		up(n);
		int diff = gh(N[n].c[0]) - gh(N[n].c[1]), d;
		if (diff >= 2) d = 0;
		else if (diff <= -2) d = 1;
		else return n;
		down(N[n].c[d]);
		if (gh(N[N[n].c[d]].c[d]) + 1 < gh(N[n].c[d]))
			N[n].c[d] = rotate(N[n].c[d], !d);
		return rotate(n, d);
	}
	static int merge_root(int l, int n, int r) {
		if (gh(l) + 1 < gh(r))
			return down(r), N[r].c[0] = merge_root(l, n, N[r].c[0]), balance(r);
		else if (gh(r) + 1 < gh(l))
			return down(l), N[l].c[1] = merge_root(N[l].c[1], n, r), balance(l);
		else return N[n].c = { l, r }, balance(n);
	}
	static std::tuple<int, int> split(int n, int k) {
		if (n != -1) down(n);
		if (k == 0) return { -1, n };
		if (k == N[n].s) return { n, -1 };
		if (k <= gs(N[n].c[0])) {
			auto [l, r] = split(N[n].c[0], k);
			return { l, merge_root(r, n, N[n].c[1]) };
		} else {
			auto [l, r] = split(N[n].c[1], k - gs(N[n].c[0]) - 1);
			return { merge_root(N[n].c[0], n, l), r };
		}
	}
	static int merge(int l, int r) {
		if (r == -1) return l;
		auto [x, nr] = split(r, 1);
		return merge_root(l, x, nr);
	}
	static void clear(int n) {
		if (n == -1) return;
		clear(N[n].c[0]);
		clear(N[n].c[1]);
		bank.push_back(n);
	}
	AVL(int v) : root(v) {}
public:
	int root;
	AVL() : root(-1) {}
	AVL& operator+=(AVL const& o) { return root = merge(root, o.root), *this; }
	friend AVL operator+(AVL a, AVL const& b) { return a += b; }
	std::tuple<AVL, AVL> split(int k) {
		auto [l, r] = split(root, k);
		return { AVL(l), AVL(r) };
	}
	void clear() { clear(root); root = -1; }
	Node& get_root() { return N[root]; }
	/*
	static AVL make_avl(Node n)
	{
		int root;
		if(bank.empty()) root = N.size(), N.push_back(std::move(n));
		else root = bank.back(), bank.pop_back(), N[root] = n;
		return AVL(root);
	} */
};
typedef AVL::Node Node;
std::vector<Node> AVL::N;
std::vector<int> AVL::bank;
