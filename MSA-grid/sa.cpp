/**
    Weighted Suffix Array
    Copyright (C) 2017 Panagiotis.Charalampopoulos, Costas.Iliopoulos, Chang.Liu, Solon.Pissis
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <map>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <boost/functional/hash.hpp>

#include "defs.h"
#include "sa.h"

#include <divsufsort.h>
#include <sdsl/rmq_support.hpp>

using namespace std;
using namespace sdsl;

void right_compacted_trie ( unordered_set<int64_t > &anchors, int64_t  * SA, int64_t  * LCP, int64_t  n, int64_t  * RSA, int64_t  * RLCP, int64_t  g )
{
	int64_t  ii = 0; //this is the index over RSA[o..g-1] and the RLCP[0..g-1], where g is the number of anchors
	int64_t  minLCP = n; //this is to maint64_t ain the minimum LCP in a range over LCP[i..j]: stores the LCP of SA[i] and SA[j]

	for( int64_t  i = 0; i < n; i++ ) // in lex order
	{

		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		auto it = anchors.find( SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		// store this suffix
			if ( ii == 0 )	RLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( SA[i-1] == RSA[ii-1] )     // if the immediately prior suffix was added
					RLCP[ii] = LCP[i];	// then the LCP value is the correct one for RLCP
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);	//otherwise, we should take the minimum in the range
			}
  			//cout<<"RSA[i]: "<< RSA[ii] <<" RLCP[i]: "<<RLCP[ii]<<"\n"; getchar();
			minLCP = n; // set this to something high to get the _FIRST_ next minimum value, because a new range STARTS
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
}

/* Constructs the left compacted trie given the anchors and the SA of the whole string in O(n) time */
void left_compacted_trie ( unordered_set<int64_t > &anchors, int64_t  * SA, int64_t  * LCP, int64_t  n, int64_t  * RSA, int64_t  * RLCP, int64_t  g )
{
	int64_t  ii = 0;
	int64_t  minLCP = n;
	for( int64_t  i = 0; i < n; i++ ) // in lex order
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		auto it = anchors.find( ( n - 1) - SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		// store this suffix
			if ( ii == 0 )	RLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( SA[i-1] == RSA[ii-1] ) // if the immediately prior suffix was added
					RLCP[ii] = LCP[i];	//then the LCP value is the correct one for RLCP
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);	//otherwise, we should take the minimum in the range
			}

			//cout<<"LSA[i]: "<< RSA[ii] <<" LLCP[i]: "<< RLCP[ii]<<"\n"; getchar();
			minLCP = n; //set this to something high to get the FIRST next minimum value
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
	if(ii < anchors.size()){
		while(ii != anchors.size()){
			RSA[ii] = n;
			ii++;
		}
	}
}

int64_t  fast_RMQ ( int64_t  * arr, int64_t  L, int64_t  R, int64_t  n, std::unordered_map<pair<int64_t ,int64_t >, int64_t , boost::hash<pair<int64_t ,int64_t > >> &rmq )
{
	if ( R - L > 1 ) // int64_t ernal nodes
	{
		int64_t  M = (L+R)/2;
		int64_t  a = fast_RMQ ( arr, L, M, n, rmq ); //Recurse on the left
		int64_t  b = fast_RMQ ( arr, M, R, n, rmq ); //Recurse on the right
		int64_t  value = std::min(a, b); 	     //This is the minimum of arr[L+1..M]

		pair<int64_t ,int64_t > p(L+1, R);
		pair<pair<int64_t ,int64_t >,int64_t > np(p, value);
		rmq.insert(np);
		return value;
	}
	else  // Base case: leaf nodes
	{
		int64_t  value;
		if ( R >= n ) value = 0;
		else value = std::min (arr[L+1], arr[R]);
		pair<int64_t ,int64_t > p(L+1, R);
		pair<pair<int64_t ,int64_t >,int64_t > np(p, value);
		rmq.insert(np);
		return value;
	}
}


int64_t  LCParray ( unsigned char * text, size_t n, int64_t  * SA, int64_t  * ISA, int64_t  * LCP )
{
        size_t i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }

        return ( 1 );
}