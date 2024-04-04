#include <cstdio>
#include <vector>

template<typename T, typename U> bool ckmax(T &a, U const& b) {return b>a?a=b,1:0;}
template<typename T, typename U> bool ckmin(T &a, U const& b) {return b<a?a=b,1:0;}

#define all(x) std::begin(x), std::end(x)
typedef long long ll;
ifdef LOCAL
template<typename T>
struct vector: std::vector<T> {
	template<typename... V>
		vector(V&&... args): std::vector<T>(std::forward<V>(args)...) { }
	T& operator[] (size_t idx) {return this->at(idx);}
	T const& operator[] (size_t idx) const {return this->at(idx);}
};
#else
using std::vector;
#endif

int main() {
}
