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

#include <iostream>
#include <string>
#include <chrono>
#include <list>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <malloc.h>

#include "defs.h"
#include "heavy_string.h"
#include "input.h"
#include "estimation.h"
#include "property_string.h"
#include "anchors_new.h"
#include "sa.h"
#include "grid.h"
#include "pattern_matching.h"

#include <divsufsort.h>
#include <sdsl/rmq_support.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

typedef grid_point point;
typedef grid_query query;

using namespace std;
using namespace std::chrono;
using get_time = chrono::steady_clock;

bool sort_sa(const pair<int,int> &a,const pair<int,int> &b)
{
       return a.first<b.first;
}

bool is_valid(vector<vector<double>>& text, string& alph, string& p, int64_t   pos, double z){
		unordered_map<char, int> mapping;
		pos = pos%(text.size());
		if(pos + p.size() > text.size()) return false;
		for(int i = 0; i < alph.size(); i++){
			mapping[alph[i]] = i;
		}
		double cum_prob = 1;
		for(int i = 0; i < p.size(); i++){
			cum_prob *= text[pos+i][mapping[p[i]]];
		}
		return (cum_prob*z >= 1)?true:false;
}

int main (int argc, char ** argv )
{
    Settings st = decode_switches(argc, argv);
	istream& text_file = st.text.is_open()?st.text:cin;
	ostream& output_file = st.output.is_open()?st.output:cout;
	ofstream result;
	
	struct mallinfo2 mi;
    mi = mallinfo2();
	double begin_ram = mi.hblkhd + mi.uordblks;
	
	double z = st.z;
	int ell = st.ell;

	string alphabet;
	vector<vector<double>> text;
	vector<vector<double>> rev_text;

	int64_t   N;
	text_file >> N;
	text_file >> alphabet;
	int64_t   A = alphabet.size();
	for (int64_t   i = 0; i < N; ++i) {
        double sum = 0;
        vector<double> symbol(A, 0);
        for (int64_t   j = 0; j < A; ++j) {
            text_file >> symbol[j];
            sum += symbol[j];
        }
        if (abs(sum-1) > EPS) {
            cerr << "Probabilities at position " << i << " do not sum up to 1" << endl;
            throw 1;
        }
        text.emplace_back(symbol);
		rev_text.emplace(rev_text.begin(), symbol);
    }
	int k = ceil(3 + log2(ell) / log2(alphabet.size()));
	int w = ell - k + 1;

	auto begin = get_time::now();
	cout << "index begin" << endl;
	Estimation fS(text, alphabet, z);
	string zstrs;
	unordered_set<int64_t> f_mini_pos;
	
	int64_t   ii = 0;
	for(PropertyString const & s : fS.strings()){
		//fT += s;
		zstrs += s.string();
		std::vector<int64_t  > M;
		minimizers_with_kr(s.string(), M,w, k);
		for(auto it : M){
			if(s._pi[it] >= k){
				f_mini_pos.emplace(it + ii*N);
			}
		}
		ii++;
	}
	
	int64_t   Nz = zstrs.size();
	int64_t   g = f_mini_pos.size();
	string rev_zstrs(zstrs.rbegin(), zstrs.rend());

	HeavyString fH(text, zstrs, alphabet);

	int64_t   * fSA		= new int64_t   [Nz];
	int64_t   * fLCP	= new int64_t   [Nz];
	int64_t   * rSA		= new int64_t   [Nz];
	int64_t   * rLCP	= new int64_t   [Nz];
	
	begin = get_time::now();
	unsigned char * seq = (unsigned char *)zstrs.c_str();
	divsufsort64( seq, fSA,  Nz );
	int64_t  * iSA = new int64_t   [Nz];
	for(int64_t  i = 0; i < Nz; i++){
		iSA[fSA[i]] = i;
	}
	LCParray( seq, Nz, fSA, iSA, fLCP );
	
	seq = (unsigned char *)rev_zstrs.c_str();
	divsufsort64( seq, rSA,  Nz );
	
	for(int64_t   i = 0; i < Nz; i++){
		iSA[rSA[i]] = i;
	}
	LCParray( seq, Nz, rSA, iSA, rLCP );
	
	fS.clear();

	int64_t   * RSA	 = new int64_t   [g];
	int64_t   * RLCP = new int64_t   [g];
	int64_t   * LSA	 = new int64_t   [g];
	int64_t   * LLCP = new int64_t   [g];
	
	right_compacted_trie ( f_mini_pos, fSA, fLCP, Nz, RSA, RLCP, g );
	left_compacted_trie ( f_mini_pos, rSA, rLCP, Nz, LSA, LLCP, g );

	delete[] fSA;
	delete[] fLCP;
	delete[] rSA;
	delete[] rLCP;
	
	vector<int> tmp_llcp(LLCP, LLCP+Nz);
	vector<int> tmp_rlcp(RLCP, RLCP+Nz);	
	rmq_succinct_sct<> lrmq ( &tmp_llcp );
	rmq_succinct_sct<> rrmq ( &tmp_rlcp );
	
	string().swap(zstrs);
	string().swap(rev_zstrs);	
	vector<int>().swap(tmp_llcp);
	vector<int>().swap(tmp_rlcp);

	vector<pair<int64_t  ,int64_t  >> l;
	vector<pair<int64_t  ,int64_t  >> r;
	for ( int64_t   i = 0; i < g; i++ )
  	{
  		l.push_back( make_pair( Nz-1-LSA[i], i) );
 		r.push_back( make_pair( RSA[i], i ) );	
	}
	
 	sort(l.begin(),l.end(),sort_sa);
	sort(r.begin(),r.end(),sort_sa);

  	std::vector<point> points;
  	  	
  	for ( int64_t   i = 0; i < g; i++ )
  	{
 
		point to_insert;
		to_insert.row = l.at(i).second+1;
		to_insert.col = r.at(i).second+1;
		
		to_insert.level = 0;
		to_insert.label = l.at(i).first;
		points.push_back(to_insert); 

  	}
  	  	
  	//constructing grid with bd-anchors
  	grid construct;
  	construct.build( points, 0 );
  	
	vector<pair<int64_t  ,int64_t  >>().swap(l);
	vector<pair<int64_t  ,int64_t  >>().swap(r);
	vector<point>().swap(points);
	
	mi = mallinfo2();
	
	double end_ram = mi.hblkhd + mi.uordblks;
	auto end = get_time::now();
	auto diff2 = end - begin;
	output_file << "Construct Time:  "<< chrono::duration_cast<chrono::milliseconds>(diff2).count()<<" ms"<<endl;	
	output_file << "Construct space:" << (end_ram-begin_ram)/1000000 << " MB" << endl;
	

	if(!st.patterns.empty()){
		string pfile = st.patterns;

		ifstream file(pfile, std::ios_base::in | std::ios_base::binary);
		boost::iostreams::filtering_istream patterns;
		patterns.push(boost::iostreams::gzip_decompressor());
		patterns.push(file);	
		begin = get_time::now();
		for (string pattern; getline(patterns, pattern); ){
			if(pattern.size() < k) continue;
	//			output_file << pattern << ":"; 
			size_t j = find_minimizer_index(pattern, k);
			string left_pattern = pattern.substr(0, j+1);
			reverse(left_pattern.begin(), left_pattern.end());
			pair<int64_t ,int64_t> left_interval = rev_pattern_matching ( left_pattern, fH, LSA, LLCP, lrmq, g ); 

			string right_pattern = pattern.substr(j);
			pair<int64_t ,int64_t> right_interval = pattern_matching ( right_pattern, fH, RSA, RLCP, rrmq, g ); 
			
			if ( left_interval.first <= left_interval.second  && right_interval.first <= right_interval.second )
			{
			
				//Finding rectangle containing bd-anchors in grid
				grid_query rectangle;
				rectangle.row1 = left_interval.first+1;
				rectangle.row2 = left_interval.second+1;
				rectangle.col1 = right_interval.first+1;
				rectangle.col2 = right_interval.second+1;			
				
				vector<size_t> result;
				construct.search_2d(rectangle, result);
			
				set<int64_t> valid_res;
				for(int64_t i = 0; i < result.size(); i++){
					if(is_valid(text, alphabet, pattern, (RSA[result.at(i)-1]-j), z)){
						valid_res.insert((RSA[result.at(i)-1]-j)%N);
					}
				}
				// for(auto i : valid_res)
					// output_file<< i << " ";
				
			}
			output_file << endl;	
		}
		end = get_time::now();
		auto diff3 = end - begin;
		output_file << pfile << " Search Time:  "<< chrono::duration_cast<chrono::milliseconds>(diff3).count()<<"ms "<<endl;
	}

	return 0;
}

