#include "../utilities/template.h"

#include "../../content/data-structures/FenwickTree2d.h"

int main() {
	rep(it,0,1000000) {
		FT2 ft(12);
		vector<tuple<int, int, int>> upd;
		int c = rand() % 20;
		rep(i,0,c) {
			upd.emplace_back(rand() % 12, rand() % 12, rand() % 10 - 5);
		}

		vector<vi> grid(12, vi(12)), sumto(13, vi(13));
		for(auto &pa: upd)
			ft.fakeUpdate(get<0>(pa), get<1>(pa));
		ft.init();
		for(auto &pa: upd) {
			grid[get<0>(pa)][get<1>(pa)] += get<2>(pa);
			ft.update(get<0>(pa), get<1>(pa), get<2>(pa));
		}

		rep(i,0,13) {
			rep(j,0,13) {
				ll v = ft.query(i, j);
				if (i == 0 || j == 0) assert(v == 0);
				else {
					sumto[i][j] = grid[i-1][j-1] + sumto[i-1][j] + sumto[i][j-1] - sumto[i-1][j-1];
					assert(v == sumto[i][j]);
				}
			}
		}
	}
	// larger / sparse coordinates, negative y, duplicate fakeUpdates, big values
	mt19937_64 rng(11);
	{ FT2 e(0); e.init(); assert(e.query(0, 5) == 0); }
	{ FT2 e(5); e.init(); rep(i,0,6) assert(e.query(i, 100) == 0); }
	rep(it,0,3000) {
		int X = (int)(rng() % 40) + 1 + (it % 100 == 0 ? 3000 : 0);
		int Y = it % 3 == 0 ? 3 : it % 3 == 1 ? 60 : 1000000000;
		ll V = it % 2 ? 10 : (ll)1e15;
		int c = (int)(rng() % 60);
		vector<tuple<int, int, ll>> upd;
		rep(i,0,c) {
			if (i && rng() % 4 == 0) upd.push_back(upd[rng() % i]), get<2>(upd.back()) = (ll)(rng() % (2 * V + 1)) - V;
			else upd.emplace_back((int)(rng() % X), (int)(rng() % (2 * (ll)Y + 1)) - Y, (ll)(rng() % (2 * V + 1)) - V);
		}
		FT2 ft(X);
		for (auto& u : upd) ft.fakeUpdate(get<0>(u), get<1>(u));
		ft.init();
		vector<tuple<int, int, ll>> done;
		for (auto& u : upd) {
			ft.update(get<0>(u), get<1>(u), get<2>(u));
			done.push_back(u);
			rep(q,0,4) {
				int x = (int)(rng() % (X + 1)), y;
				if (rng() % 2) y = get<1>(upd[rng() % c]) + (int)(rng() % 3) - 1;
				else y = (int)(rng() % (2 * (ll)Y + 3)) - Y - 1;
				if (rng() % 8 == 0) x = X;
				if (rng() % 8 == 0) y = INT_MAX;
				if (rng() % 16 == 0) y = INT_MIN;
				ll sum = 0;
				for (auto& d : done) if (get<0>(d) < x && get<1>(d) < y) sum += get<2>(d);
				assert(ft.query(x, y) == sum);
			}
		}
	}
	cout<<"Tests passed!"<<endl;
}
