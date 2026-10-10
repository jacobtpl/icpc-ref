/**
 * Author: Siyong
 * Date: 2023
 * Description: FWHT (in place). The size of $a$ must be a power of two ($\ge 1$).
 * For inverse fwht: call fwht and divide all elements by (2 ** L) (size of array)
 * Time: O(N \log N)
 */
template<typename T>
void fwht(vector<T>& a)
{
	int L=__builtin_ctz(a.size()); // should be equal to 31 - __builtin_clz(a.size());
	assert(L == 31 - __builtin_clz(a.size()));
 
	for(int i=0;i<L;++i)
		for(int j=0;j<1<<L;j+=1<<i+1)
			for(int k=0;k<1<<i;++k)
			{
				T v=a[j|1<<i|k];
				a[j|1<<i|k] = a[j|k]-v;
				a[j|k] += v;
			}
}
