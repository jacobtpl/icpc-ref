/**
 * Author: Benq
 * Description: trie of right-half of all palindromic substrings
 * link points to longest suffix-palindrome of node
 * oc stores num occurences after resolveOc()
 * ans[i][b] is min factorization of prefix i (1-indexed) with parity b
 * Time: O(n\log n)
 * Status: tested
 */
struct PalTree {
    static const int ASZ = 26;
    struct node {
        array<int,ASZ> to = array<int,ASZ>();
        int len, link, oc = 0; // # occurrences of pal
        int slink = 0, diff = 0;
        array<int,2> seriesAns;
        node(int _len, int _link) : len(_len), link(_link) {}
    };
    string s = "@"; vector<array<int,2>> ans = {{0,INF}};
    vector<node> d = {{0,1},{-1,0}}; // dummy pals of len 0,-1
    int last = 1;
    int getLink(int v) {
        while (s[sz(s)-d[v].len-2] != s.back()) v = d[v].link;
        return v;
    }
    void updAns() { // serial path has O(log n) vertices
        ans.pb({INF,INF});
        for (int v = last; d[v].len > 0; v = d[v].slink) {
            d[v].seriesAns=ans[sz(s)-1-d[d[v].slink].len-d[v].diff];
            if (d[v].diff == d[d[v].link].diff) 
                rep(i,0,2) ckmin(d[v].seriesAns[i], d[d[v].link].seriesAns[i]);
            // start of previous oc of link[v]=start of last oc of v
            rep(i,0,2) ckmin(ans.back()[i],d[v].seriesAns[i^1]+1);
        }
    }
    void addChar(char C) {
        s += C; int c = C-'a'; last = getLink(last);
        if (!d[last].to[c]) {
            d.emplace_back(d[last].len+2,d[getLink(d[last].link)].to[c]);
            d[last].to[c] = sz(d)-1;
            auto& z = d.back(); z.diff = z.len-d[z.link].len;
            z.slink = z.diff == d[z.link].diff ? d[z.link].slink : z.link;
        } // max suf with different dif
        last = d[last].to[c]; ++d[last].oc;
        updAns();
    }
    void resolveOc() { for (int i=sz(d)-1;i>=2;i--) d[d[i].link].oc += d[i].oc; }
};