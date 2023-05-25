/**
 * Author: Siyong Huang
 * Date: 2023-05-01
 * License: CC0
 * Source: My head
 * Description: what it says
 * Usage: just do it
 * Status: prob works
 */
#include <vector>

template<int B = 20>
struct BigUInt: public std::vector<int>
{
public:
	static int const ALL = (1 << B) - 1;
	BigUInt(ll x = 0): vector()
	{
		for(;x>0;x>>=B)
			emplace_back(x & ALL);
		check();
	}

	void check()
	{
#ifdef LOCAL
		assert(empty() || back() > 0);
		for(int x: *this) assert(0 <= x);
		for(int x: *this) assert(x < (1 << B));
#endif
	}
	BigUInt& operator+= (BigUInt const& o)
	{
		int carry = 0;
		for(int i=0;carry || i < o.size();++i)
		{
			if(i < o.size()) carry += o[i];
			if(i < size()) at(i) += carry;
			else emplace_back(carry & ALL);
			carry = at(i) >> B;
			at(i) &= ALL;
		}
		check();
		return *this;
	}

	BigUInt& operator-= (BigUInt const& o)
	{
		int carry = 0;
		for(int i=0;carry || i < o.size();++i)
		{
			if(i < o.size()) carry += o[i];
			if(i < size()) at(i) -= carry;
			else emplace_back(-carry);
			assert(back() >= 0);
			if(at(i) < 0) carry = 1, at(i) += (1 << B); else carry = 0;
			assert(0 <= at(i) && at(i) < (1 << B));
		}
		for(;!empty() && back() == 0;) pop_back();
		check();
		return *this;
	}

	friend BigUInt operator+ (BigUInt a, BigUInt const& b) {return a+=b;}
	friend BigUInt operator- (BigUInt a, BigUInt const& b) {return a-=b;}

	bool operator< (BigUInt const& o) {
		if(size() != o.size()) return size() < o.size();
		for(int i=(int)size()-1;i>=0;--i)
			if(at(i) != o[i])
				return at(i) < o[i];
		return 0;
	}
	bool operator> (BigUInt const& o) {
		if(size() != o.size()) return size() > o.size();
		for(int i=(int)size()-1;i>=0;--i)
			if(at(i) != o[i])
				return at(i) > o[i];
		return 0;
	}

	BigUInt operator<< (uint shift)
	{
		uint full = shift / B;
		uint partial = shift - full * B;
		BigUInt o;
		o.assign(full + 1, 0);
		for(int i=0;i<size();++i)
		{
			o[i + full] |= at(i) << partial & ALL;
			o.push_back(at(i) >> B - partial);
		}
		for(;!o.empty() && o.back() == 0;o.pop_back());
		o.check();
		return o;
	}
	BigUInt operator>> (uint shift)
	{
		uint full = shift / B;
		uint partial = shift - full * B;
		BigUInt o;
		for(int i=full;i<size();++i)
		{
			if(i-full >= 0)
				o.v[i - full - 1] |= at(i) << (B - partial) & ALL;
			o.push_back(at(i) >> partial);
		}
		for(;!o.empty() && o.back() == 0;) o.pop_back();
		o.check();
		return o;
	}

	BigUInt& operator>>= (uint shift) { return *this = *this >> shift; }
	BigUInt& operator<<= (uint shift) { return *this = *this << shift; }
};
