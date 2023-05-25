/**
 * Author: Siyong Huang
 * Date: 2023-05-01
 * License: CC0
 * Source: My head
 * Description: what it says
 * Usage: just do it
 * Status: not super optimized
 */
#include <vector>
#include <map>
struct SA
{
public:
	std::vector<std::map<char, int> > adj;
	std::vector<int> link, dis;
	int N;
	SA(): adj(1), link(1, -1), dis(1, 0), N(1) {}

	int new_node(int v=-1)
	{
		if(v == -1)
			adj.emplace_back(), link.emplace_back(), dis.emplace_back();
		else
			adj.push_back(adj[v]), link.push_back(link[v]), dis.push_back(dis[v]);
		return N++;
	}

	int go(int p, int c)
	{
		auto it = adj[p].find(c);
		if(dis[it->second] == dis[p] + 1)
			return it->second;
		else
		{
			int q = it->second, n = new_node(q);
			dis[n] = dis[p] + 1, link[q] = n;
			for(;p != -1 && (it = adj[p].find(c))->second == q;p = link[p])
				it->second = n;
			return n;
		}
	}

	int append(int p, char c)
	{
		auto it = adj[p].find(c);
		if(it != adj[p].end())
			return go(p, c);

		int n = new_node();
		dis[n] = dis[p] + 1;
		for(;p != -1 && adj[p].find(c) == adj[p].end();p = link[p])
			adj[p].insert({c, n});
		if(p == -1)
			link[n] = 0;
		else
			link[n] = go(p, c);
		return n;
	}

	int add(std::string const &s)
	{
		int n = 0;
		for(char c: s)
			n = append(n, c);
		return n;
	}
};
