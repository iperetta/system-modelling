/* code by: iperetta@ieee.org */
#ifndef CPRNGXS_H_
#define CPRNGXS_H_

#include <cmath>
#include <vector>
/* Adapted from http://xorshift.di.unimi.it  *
 * (CC0 1.0) 2014, by Sebastiano Vigna       */

#ifndef UINT64_MAX
#define UINT64_MAX 0xffffffffffffffff
typedef unsigned long long uint64_t;
#endif

class cPRNGXS_64 { // xorshift64*
  private:
    uint64_t x;
  public:
    cPRNGXS_64(uint64_t seed) : x(seed)
    {
		//// MurmurHash3's 64-bit finalizer =======
		//x ^= x >> 33;
		//x *= 0xff51afd7ed558ccd;
		//x ^= x >> 33;
		//x *= 0xc4ceb9fe1a85ec53;
		//x ^= x >> 33;
		//// ======================================
	};
	inline uint64_t next() 
	{
		x ^= x >> 12; // a
		x ^= x << 25; // b
		x ^= x >> 27; // c
		return x * 2685821657736338717LL;
	}
    ~cPRNGXS_64() {};    
};

class cPRNGXS_1024 { // xorshift1024*
  private:
    uint64_t s[ 16 ]; 
	int p;
	bool ready;
	double next_normal;
  public:
    cPRNGXS_1024(uint64_t seed) : p(0), ready(false), next_normal(0.)
    {
		cPRNGXS_64 randomseed(seed);
		for(int i = 0; i < 16; i++)
			s[i] = randomseed.next();
	};
	inline uint64_t next() 
	{
		uint64_t s0 = s[ p ];
		uint64_t s1 = s[ p = ( p + 1 ) & 15 ];
		s1 ^= s1 << 31; // a
		s[ p ] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
		return s[ p ] * 1181783497276652981LL;
	}
	int integer(const int &a = 0, const int &b = 100)
	{
		return a + round((b - a)*static_cast<double>(next())
			/static_cast<double>(UINT64_MAX));
	}
	double real(const double &a = 0., const double &b = 1.)
	{
		return a + (b - a)*static_cast<double>(next())
			/static_cast<double>(UINT64_MAX);
	}
	double normal(const double &mu = 0., const double &std = 1.)
	{
		if(ready)
		{
			ready = false;
			return mu + std*next_normal;
		}
		else
		{
			// Box-Muller transformation
			double u1, u2, rsq, fac;
			do
			{
				u1 = -1. + 2.*static_cast<double>(next())
					/static_cast<double>(UINT64_MAX);
				u2 = -1. + 2.*static_cast<double>(next())
					/static_cast<double>(UINT64_MAX);
				rsq = u1*u1 + u2*u2;
			}
			while(rsq >= 1. || rsq == 0.);
			fac = sqrt(-2.*log(rsq)/rsq);
			ready = true;
			next_normal = u2*fac;
			return mu + std*u1*fac;
		}
	}
	std::vector<size_t> randperm(const size_t &N)
	{
		// return list with 0..N-1 randomly permutated
		size_t idx1, idx2, aux;
		std::vector<size_t> rp;
		for(size_t i = 0; i < N; i++)
			rp.push_back(i);
		for(size_t i = 0; i < N; i++)
		{
			idx1 = (size_t)this->integer(0,N-1);
			idx2 = (size_t)this->integer(0,N-1);
			aux = rp.at(idx1);
			rp.at(idx1) = rp.at(idx2);
			rp.at(idx2) = aux;
		}
		return rp;
	}
    ~cPRNGXS_1024() {};    
};

#endif
