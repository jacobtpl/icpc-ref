#include "../utilities/template.h"
#include "../utilities/utils.h"
#include "../utilities/random.h"

#include "../../content/graph/WeightedMatching.h"
#include <bits/extc++.h> /// include-line, keep-include
template<class T> bool ckmin(T& a, const T& b) {
    if (a <= b) return false;
    a = b;
    return true;
}
#include "../../content/graph/MinCostMaxFlow.h"

void test(int N, int mxCost, int iters) {
    for (int it = 0; it < iters; it++) {
        int n = randIncl(0, N), m = randIncl(0, N);
        if (n > m)
            swap(n, m);

        MCMF mcmf(n + m + 2);
        int s = 0;
        int t = 1;
        for (int i = 0; i < n; i++)
            mcmf.addEdge(s, i + 2, 1, 0);
        for (int i = 0; i < m; i++)
            mcmf.addEdge(2 + n + i, t, 1, 0);

        vector<vi> cost(n, vi(m));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                cost[i][j] = randRange(-mxCost, mxCost);
                mcmf.addEdge(i + 2, 2 + n + j, 1, cost[i][j]);
            }
        }
        mcmf.setpi(s);
        auto maxflow = mcmf.maxflow(s, t);
        auto matching = hungarian(cost);
        assert(maxflow.first == n);
        assert(maxflow.second == matching.first);
        int matchSum = 0;
        set<int> used;
        for (int i = 0; i < n; i++) {
            matchSum += cost[i][matching.second[i]];
            assert(used.count(matching.second[i]) == 0);
            used.insert(matching.second[i]);
        }
        assert(matchSum == matching.first);
    }
}
template<class T>
void verifyBrute(const vector<vector<T>>& cost) {
    int n = sz(cost), m = n ? sz(cost[0]) : 0;
    vi cols(m);
    iota(all(cols), 0);
    T best = numeric_limits<T>::max();
    do {
        T sum = 0;
        rep(i,0,n) sum += cost[i][cols[i]];
        best = min(best, sum);
    } while (next_permutation(all(cols)));
    auto [value, matching] = hungarian(cost);
    assert(value == best && sz(matching) == n);
    T sum = 0;
    set<int> used;
    rep(i,0,n) {
        assert(0 <= matching[i] && matching[i] < m);
        assert(used.insert(matching[i]).second);
        sum += cost[i][matching[i]];
    }
    assert(sum == value);
}

signed main() {
    assert(hungarian({}).first == 0);
    assert(hungarian({{3, 1}, {2, 4}}).first == 3);
    const vector<vector<ll>> large = {{4000000000LL, -5000000000LL},
                                     {-6000000000LL, 7000000000LL}};
    verifyBrute(large);
    mt19937_64 rng(17);
    rep(it,0,2000) {
        int n = int(rng() % 5), m = n + int(rng() % (7 - n));
        vector<vector<ll>> cost(n, vector<ll>(m));
        vector<vector<double>> realCost(n, vector<double>(m));
        rep(i,0,n) rep(j,0,m) {
            cost[i][j] = ll(rng() % 2000000000000001ULL) - 1000000000000000LL;
            realCost[i][j] = double(int(rng() % 81) - 40) / 4;
        }
        verifyBrute(cost);
        verifyBrute(realCost);
    }
    test(25, 5, 1000);
    test(100, 1000, 100);
    test(100, 1, 50);
    test(5, 5, 10000);
    cout << "Tests passed!" << endl;
}
