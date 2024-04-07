/**
 * Author: Siyong
 * Date: 2023
 * License: CC0
 * Description: Prime sieve but slow, for generating all primes smaller than LIM.
 * Time: LIM=1e9 $\approx$ 8.5s
 * Status: kinda tested
 */
int const LIM = 1e7+5;
vector<bool> cp;
vi pr, nx, lp, cnt;
void sieve()
{
	cp.assign(LIM, 0), nx.assign(LIM, -1), lp.assign(LIM, -1), cnt.assign(LIM, -1);
	for(int i=2;i<LIM;++i) {
		if(!cp[i])
			lp[i] = pr.size(), nx[i] = cnt[i] = 1, pr.push_back(i);
		for(int j=0,k;j<pr.size() && (k=i*pr[j])<LIM;++j) { // pr[j]<(LIM+i-1)/i, if there's overflow
			cp[k] = 1, lp[k] = j;
			if(j == lp[i]) {
				nx[k] = nx[i], cnt[k] = cnt[i]+1; break;
			} else nx[k] = i, cnt[k] = 1;;
		}
	}
}
/*
int main() {
	sieve();
	int N; scanf("%d", &N);
	for(int i=0;i<N;++i) {
		int x;
		scanf("%d", &x);
		for(;x>1;x=nx[x])
			for(int i=0;i<cnt[x];++i)
				printf("%d ", pr[lp[x]]);
		printf("\n");
	}
	return 0;
}*/
