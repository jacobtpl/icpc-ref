#include "../utilities/template.h"

#include "../../content/strings/Zfunc.h"

template <class F> void gen(string &s, int at, int alpha, F f) {
    if (at == sz(s))
        f();
    else {
        rep(i, 0, alpha) {
            s[at] = (char)('a' + i);
            gen(s, at + 1, alpha, f);
        }
    }
}

void test(const string &s) {
    int n = sz(s);
    vi found = Z(s);
    vi expected(n, 0);
    rep(i, 1, n) { // exclude index 0 (!)
        int j = 0;
        while (i + j < n && s[i + j] == s[j])
            j++;
        expected[i] = j;
    }
    assert(found == expected);
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    rep(n, 0, 13) {
        string s(n, 'x');
        gen(s, 0, 3, [&]() { test(s); });
    }
    rep(n, 0, 11) {
        string s(n, 'x');
        gen(s, 0, 4, [&]() { test(s); });
    }
    // longer random / periodic strings, arbitrary bytes
    mt19937 rng(555);
    rep(it, 0, 100000) {
        int n = (int)(rng() % 80), kind = (int)(rng() % 3);
        int alpha = kind == 0 ? 256 : 1 + (int)(rng() % 3);
        int p = 1 + (int)(rng() % 6);
        string s(n, 0);
        rep(i, 0, n) s[i] = (char)(rng() % (unsigned)alpha);
        if (kind == 2) rep(i, p, n) if (rng() % 30) s[i] = s[i - p];
        test(s);
    }
    cout<<"Tests passed!"<<endl;
}
