/**
 * Author: jacobtpl
 * Date: 2024
 * Source: Benq
 * Description: Aho-Corasick with unfixed, trie of all strings
 * Link stores node of longest suffix that is also in the trie
 * Only works after adding all strings and calling build().
 * Status: tested
 */
struct AhoCorasick {
    struct node {
        unordered_map<char,int> next; // can use array
        int link = -1;
        int output = 0; // adapt as necessary
        node() {}
    };
    vector<node> nodes{1};
    void add_string(const string& s) {
        int v = 0;
        for (char c:s) {
            if (nodes[v].next.find(c) == nodes[v].next.end()) {
                nodes[v].next[c] = sz(nodes);
                nodes.emplace_back();
            }
            v = nodes[v].next[c];
        }
        nodes[v].output = 1; // update output
    }
    int nx(int j, char c) {
        while (j != -1 && !nodes[j].next.count(c)) j = nodes[j].link;
        return (j == -1) ? 0 : nodes[j].next[c];
    }
    void build() {
        queue<int> q; q.push(0);
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (auto p : nodes[v].next) {
                int x = p.second;
                nodes[x].link = nx(nodes[v].link, p.first);
                // update output if needed
                q.push(x);
            }
        }
    }
    int match(const string &text) {
        int n = sz(text), v = 0;
        for (int i = 0; i < n; i++) {
            v = nx(v, text[i]);
            if (nodes[v].output) return i;
        }
        return -1;
    }
};