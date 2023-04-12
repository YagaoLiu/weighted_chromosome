#include <unordered_set>
#include <math.h>
#include <string>

#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>

#include "pattern_matching.h"

#include <boost/functional/hash.hpp>
using namespace boost;
using namespace std;
using namespace sdsl;

int64_t  find_minimizer_index(string s, int64_t  k) {
  int64_t  minimizer_index = 0;
  for (int64_t  i = 1; i <= s.length() - k; i++) {
    if (s.substr(i, k) < s.substr(minimizer_index, k)) {
      minimizer_index = i;
    }
  }

  return minimizer_index;
}

/* Computes the length of lcp of two suffixes of two strings */
int64_t  lcp ( HeavyString & x, int64_t  M, string & y, int64_t  l )
{
	int64_t  xx = x.size();
	if ( M >= xx ) return 0;
	int64_t  yy = y.size();
	if ( l >= yy ) return 0;

	int64_t  i = 0;
	while ( ( M + i < xx ) && ( l + i < yy ) )
	{
		if ( x[M+i] != y[l+i] )	break;
		i++;
	}
	return i;
}

/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<int64_t ,int64_t > pattern_matching ( string & w, HeavyString & a, int  * SA, int  * LCP, rmq_succinct_sct<> &rmq, int64_t  n )
{
	int64_t  m = w.size(); //length of pattern
	int64_t  N = a.size(); //length of string 
	int64_t  d = -1;
	int64_t  ld = 0;
	int64_t  f = n;
	int64_t  lf = 0;
	
	pair<int64_t ,int64_t > range;
	
	while ( d + 1 < f )
	{
		int64_t  i = (d + f)/2;

		/* lcp(i,f) */
		int64_t  lcpif;
		lcpif = LCP[rmq(i+1, f)];
		
		/* lcp(d,i) */
		int64_t  lcpdi; 
		lcpdi = LCP[rmq(d+1, i)];

		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}	
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			int64_t  l = std::max (ld, lf);
			l = l + lcp ( a, SA[i] + l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				int64_t  e = i;
				while ( d + 1 < e )
				{
					int64_t  j = (d + e)/2;
					
					/* lcp(j,e) */
					int64_t  lcpje;
					lcpje = LCP[rmq(j+1, e)];

					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				int64_t  lcpde;
				lcpde = LCP[rmq(d+1, e)];

				if ( lcpde >= m )	d = std::max (d-1,( int64_t  ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					int64_t  j = (e + f)/2;
					
					/* lcp(e,j) */
					int64_t  lcpej;
					lcpej = LCP[rmq(e+1, j)];

					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				int64_t  lcpef;
				lcpef = LCP[rmq(e+1, f)];

				if ( lcpef >= m )	f = std::min (f+1,n);

				range.first = d + 1;
				range.second = f - 1;
				return range;
				
				
			}
			else if ( ( l == N - SA[i] ) || ( ( SA[i] + l < N ) && ( l != m ) && ( a[SA[i]+l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}	

		}
	}

	range.first = d + 1;
	range.second = f - 1;
	return range;
}

/* Computes the length of lcs of two suffixes of two strings */
int64_t  lcs ( HeavyString & x, int64_t  M, string & y, int64_t  l )
{
	if ( M < 0 ) return 0;
	int64_t  yy = y.size();
	if ( l >= yy ) return 0;

	int64_t  i = 0;
	while ( ( M - i >= 0 ) && ( l + i < yy ) )
	{
		if ( x[M-i] != y[l+i] )	break;
		i++;
	}
	return i;
}


/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<int64_t ,int64_t > rev_pattern_matching ( string & w, HeavyString & a, int  * SA, int  * LCP, rmq_succinct_sct<> &rmq, int64_t  n )
{
	int64_t  m = w.size(); //length of pattern
	int64_t  N = a.size(); //length of string 
	int64_t  d = -1;
	int64_t  ld = 0;
	int64_t  f = n;
	int64_t  lf = 0;
	
	pair<int64_t ,int64_t > range;
	
	while ( d + 1 < f )
	{
		int64_t  i = (d + f)/2;
		int64_t  revSA = N - 1 - SA[i];

		/* lcp(i,f) */
		int64_t  lcpif;
		lcpif = LCP[rmq(i+1, f)];
		
		/* lcp(d,i) */
		int64_t  lcpdi;
		lcpdi = LCP[rmq(d+1, i)];

		//std::cout << revSA << " " << i << " " << d << " " << f << " " << ld << " " << lf << " " << lcpdi << " " << lcpif << std::endl;
		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) {
			f = i;
		}
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}	
		else if ( ( lf <= ld ) && ( ld < lcpdi ) ){
			d = i;
		}
		else
		{
			int64_t  l = std::max (ld, lf);
			l = l + lcs ( a, revSA - l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				int64_t  e = i;
				while ( d + 1 < e )
				{
					int64_t  j = (d + e)/2;
					
					/* lcp(j,e) */
					int64_t  lcpje;
					lcpje = LCP[rmq(j+1, e)];

					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				int64_t  lcpde;
				lcpde = rmq(d+1, e);

				if ( lcpde >= m )	d = std::max (d-1,( int64_t  ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					int64_t  j = (e + f)/2;
					
					/* lcp(e,j) */
					int64_t  lcpej;
					lcpej = LCP[rmq(e+1, j)];

					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				int64_t  lcpef;
				lcpef = LCP[rmq(e+1, f)];

				if ( lcpef >= m )	f = std::min (f+1,n);

				range.first = d + 1;
				range.second = f - 1;
				return range;
				
				
			}
			else if ( ( l == N - SA[i] ) || ( ( revSA - l >= 0 ) && ( l != m ) && ( a[revSA - l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}	

		}
	}

	range.first = d + 1;
	range.second = f - 1;
	return range;
}