#include "../utilities/template.h"

#include "../../content/various/FastKnapsack.h"

int naive(vi w, int t) {
    vector<bool> can_reach(t+1);
    can_reach[0] = true;
    for (int x : w) {
        for (int i = t-x; i >= 0; --i) {
            if (can_reach[i]) can_reach[i+x] = true;
        }
    }
    for (int i = t;; i--)
        if (can_reach[i]) return i;
    assert(false);
}

int naiveBitset(const vi& w, int t) {
    static bitset<60003> bs;
    bs.reset(); bs[0] = 1;
    for (int x : w) bs |= bs << x;
    for (int i = t;; i--)
        if (bs[i]) return i;
}

// Compile with -DBENCH for timings.
int main() {
    // Edge cases
    assert(knapsack({}, 0) == 0);
    assert(knapsack({}, 100) == 0);
    assert(knapsack({0}, 0) == 0);
    assert(knapsack({0, 0, 0}, 5) == 0);
    assert(knapsack({5}, 0) == 0);
    assert(knapsack({5}, 4) == 0);
    assert(knapsack({5}, 5) == 5);
    assert(knapsack({5}, 6) == 5);
    assert(knapsack({0, 5, 0}, 4) == 0);
    assert(knapsack({3, 3, 3}, 100) == 9);
    assert(knapsack({3, 3, 3}, 8) == 6);
    assert(knapsack({1000000, 1000000}, 1999999) == 1000000);
    assert(knapsack({1000000, 1000000}, 2000000) == 2000000);
    assert(knapsack({7, 1000000, 3}, 999999) == 10);
    assert(knapsack(vi(100000, 1), 77777) == 77777);
    assert(knapsack(vi(1000, 1000), 999999) == 999000);

    const int MAX_N = 10;
    const int MAX_W = 10;
    const int iters = 1000000;
    rep(it,0,iters) {
        int n = rand() % MAX_N;
        int maxw = rand() % MAX_W + 1;
        vi w(n);
        rep(i,0,n)
            w[i] = rand()%(maxw+1);
        int t = rand() % (MAX_N*maxw);
        assert(naive(w,t) == knapsack(w,t));
    }
    // Larger cases: few distinct weights, sorted orders, targets around/above the total sum.
    rep(it,0,20000) {
        int n = rand() % 60 + 1;
        int maxw = rand() % 3 ? rand() % 30 + 1 : rand() % 1000 + 1;
        int distinct = rand() % 2 ? 1000000 : rand() % 3 + 1;
        vi pool(distinct < 10 ? distinct : 0), w(n);
        for (int& x : pool) x = rand() % (maxw + 1);
        rep(i,0,n) w[i] = pool.empty() ? rand() % (maxw + 1) : pool[rand() % distinct];
        int ord = rand() % 4;
        if (ord == 1) sort(all(w));
        if (ord == 2) sort(all(w), greater<int>());
        int sum = accumulate(all(w), 0), t;
        int kind = rand() % 4;
        if (kind == 0) t = rand() % (sum + 1);
        else if (kind == 1) t = sum + rand() % 3;
        else if (kind == 2) t = max(0, sum / 2 + rand() % 5 - 2);
        else t = rand() % (maxw + 1);
        assert(naiveBitset(w, t) == knapsack(w, t));
    }

#ifdef BENCH
    {
        auto bench = [&](const char* name, vi w, int t) {
            auto t0 = chrono::steady_clock::now();
            int r = knapsack(w, t);
            auto t1 = chrono::steady_clock::now();
            cerr << name << " N=" << sz(w) << " maxw=" << *max_element(all(w)) << " t=" << t << ": "
                << chrono::duration<double>(t1 - t0).count() << " s (answer " << r << ")" << endl;
        };
        for (pii p : {pii(100000, 1000), pii(10000, 10000), pii(1000, 100000), pii(1000000, 100)}) {
            int n = p.first, m = p.second;
            vi w(n);
            ll sum = 0;
            rep(i,0,n) sum += w[i] = rand() % m + 1;
            w[0] = m;
            ll wsum = sum;
            bench("random, t=sum/2", w, (int)(sum / 2));
            bench("random, t=sum-1", w, (int)(sum - 1));
            bench("random, t=maxw", w, m);
            // even weights with an odd target: the answer t is never reached
            vi e(n);
            sum = 0;
            rep(i,0,n) sum += e[i] = (rand() % (m / 2) + 1) * 2;
            e[0] = m;
            bench("all even, odd t", e, (int)(sum / 4 * 2 + 1));
            vi big(n, m); rep(i,0,n) if (i % 2) big[i] = m - 1;
            bench("weights m and m-1", big, (int)((ll)n * m / 2 - m / 2));
            sort(all(w));
            bench("sorted ascending", w, (int)(wsum / 2));
            reverse(all(w));
            bench("sorted descending", w, (int)(wsum / 2));
        }
    }
#endif
    cout<<"Tests passed!"<<endl;
}
