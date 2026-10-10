#include <cstdio>
#include <vector>

template<typename T, typename U> bool ckmax(T &a, U const& b) {return b>a?a=b,1:0;}
template<typename T, typename U> bool ckmin(T &a, U const& b) {return b<a?a=b,1:0;}

#define all(x) std::begin(x), std::end(x)
typedef long long ll;
#ifdef LOCAL
template<typename T>
struct vector: std::vector<T> {
	using std::vector<T>::vector;
	vector(const std::vector<T>& v) : std::vector<T>(v) {}
	vector(std::vector<T>&& v) : std::vector<T>(std::move(v)) {}
	decltype(auto) operator[] (size_t idx) {return this->at(idx);}
	decltype(auto) operator[] (size_t idx) const {return this->at(idx);}
};
#else
using std::vector;
#endif

int main() {
}
