struct PAVL {
public:
	struct Node {
	public:
		int t;
		int s, h; /* Customize */
		std::array<int, 2> c;
		int val;
		Node() : t(), s(1), h(1), c{ -1, -1 } { t = T; }
		Node(int _val) : Node() { val = _val; }
		void up() {}; /* Customize */
		void down() {};
	};
private:
	static std::vector<Node> N;
	static int T;
	static int clone(int n) {
		if (n == -1) return -1; //assert(N[n].t >= t);
		if (N[n].t == T) return n;
		return N.push_back(N[n]), N.back().t = T, N.size() - 1;
	}
	static int gh(int n) { return n != -1 ? N[n].h : 0; }
	static int gs(int n) { return n != -1 ? N[n].s : 0; }
	static void up(int n) {
		N[n].h = std::max(gh(N[n].c[0]), gh(N[n].c[1])) + 1;
		N[n].s = gs(N[n].c[0]) + gs(N[n].c[1]) + 1;
		N[n].up();
	}
	static int down(int n) { n = clone(n); return N[n].down(), n; }
	static int rotate(int n, int d) {
		n = clone(n); int o = down(N[n].c[d]);
		N[n].c[d] = N[o].c[!d], N[o].c[!d] = n;
		up(n), up(o);
		return o;
	}
	static int balance(int n) {
		assert(N[n].t == T); up(n);
		int diff = gh(N[n].c[0]) - gh(N[n].c[1]), d;
		if (diff >= 2) d = 0;
		else if (diff <= -2) d = 1;
		else return n;
		N[n].c[d] = down(N[n].c[d]);
		if (gh(N[N[n].c[d]].c[d]) + 1 < gh(N[n].c[d]))
			N[n].c[d] = rotate(N[n].c[d], !d);
		return rotate(n, d);
	}
	static int merge_root(int l, int n, int r) {
		if (gh(l) + 1 < gh(r))
			return r = down(r), N[r].c[0] = merge_root(l, n, N[r].c[0]), balance(r);
		else if (gh(r) + 1 < gh(l))
			return l = down(l), N[l].c[1] = merge_root(N[l].c[1], n, r), balance(l);
		else return N[n].c = { l, r }, balance(n);
	}
	static std::tuple<int, int> split(int n, int k) {
		if (n != -1) n = down(n);
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
		return merge_root(l, clone(x), nr);
	}
	PAVL(int v) : root(v) {}
public:
	int root;
	PAVL() : root(-1) {}
	PAVL(Node&& n) : root(N.size()) { N.push_back(n); }
	friend PAVL operator+(PAVL a, PAVL b) { return merge(a.root, b.root); }
	std::tuple<PAVL, PAVL> split(int k) {
		auto [l, r] = split(root, k);
		return { PAVL(l), PAVL(r) };
	}
	PAVL step() { ++T; return clone(root); }
	Node& get_root() { return N[root]; }
};
typedef PAVL::Node Node;
std::vector<Node> PAVL::N;
int PAVL::T;
