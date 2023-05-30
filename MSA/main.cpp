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

// #include "defs.h"
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

	int   N;
	text_file >> N;
	text_file >> alphabet;
	int   A = alphabet.size();
	for (int   i = 0; i < N; ++i) {
        double sum = 0;
        vector<double> symbol(A, 0);
        for (int   j = 0; j < A; ++j) {
            text_file >> symbol[j];
            sum += symbol[j];
        }
        if (abs(sum-1) > EPS) {
            cerr << "Probabilities at position " << i << " do not sum up to 1" << endl;
            throw 1;
        }
        text.emplace_back(symbol);
    }
	
	int k = ceil(4 * log2(ell) / log2(alphabet.size()));
	int w = ell - k + 1;

	auto begin = get_time::now();
	cout << "index begin" << endl;
	Estimation fS(text, alphabet, z);
	string zstrs;
	vector<int> mini_pos;
	
	int  ii = 0;
	for(PropertyString const & s : fS.strings()){
		zstrs += s.string();
		std::vector<int> M;
		minimizers_with_kr(s.string(), M,w, k);
		for(auto it : M)
			if(s._pi[it] >= k){
				mini_pos.push_back(it + ii*N);
			}
		ii++;
	}

	int   Nz = zstrs.size();
	string rev_zstrs(zstrs.rbegin(), zstrs.rend());

	vector<int> le;
	vector<int> re;
	extention(text, zstrs, alphabet, le, re, z);
	
	mini_pos.erase(remove_if(mini_pos.begin(), mini_pos.end(), [&](int i) { return re[i] < k; }), mini_pos.end());
	unordered_set<int> f_mini_pos(mini_pos.begin(), mini_pos.end());
	int g = f_mini_pos.size();

	HeavyString fH(text, zstrs, alphabet, f_mini_pos, le, re, true);
	
	fS.clear();
	vector<vector<double>>().swap(text);


	int   * fSA		= new int   [Nz];
	int   * fLCP	= new int   [Nz];
	int   * rSA		= new int   [Nz];
	int   * rLCP	= new int   [Nz];
	
	unsigned char * seq = (unsigned char *)zstrs.c_str();
	divsufsort( seq, fSA,  Nz );
	int  * iSA = new int   [Nz];
	for(int  i = 0; i < Nz; i++){
		iSA[fSA[i]] = i;
	}
	LCParray( seq, Nz, fSA, iSA, fLCP );

	seq = (unsigned char *)rev_zstrs.c_str();
	divsufsort( seq, rSA,  Nz );
	
	for(int   i = 0; i < Nz; i++){
		iSA[rSA[i]] = i;
	}
	LCParray( seq, Nz, rSA, iSA, rLCP );
	delete[] iSA;

	int   * RSA	 = new int   [g];
	int   * RLCP = new int   [g];
	int   * LSA	 = new int   [g];
	int   * LLCP = new int   [g];
	
	right_compacted_trie ( f_mini_pos, fSA, fLCP, Nz, RSA, RLCP, g );
	left_compacted_trie ( f_mini_pos, rSA, rLCP, Nz, LSA, LLCP, g );

	vector<int> le_r(le.rbegin(), le.rend());	
	union_find_resort ( RSA, RLCP, re, g );
	union_find_resort ( LSA, LLCP, le_r, g );

	delete[] fSA;
	delete[] fLCP;
	delete[] rSA;
	delete[] rLCP;
	
	vector<int> tmp_llcp(LLCP, LLCP+g);
	vector<int> tmp_rlcp(RLCP, RLCP+g);	
	rmq_succinct_sct<> lrmq ( &tmp_llcp );
	rmq_succinct_sct<> rrmq ( &tmp_rlcp );
	
	string().swap(zstrs);
	string().swap(rev_zstrs);	
	vector<int>().swap(tmp_llcp);
	vector<int>().swap(tmp_rlcp);
	vector<int>().swap(mini_pos);
	unordered_set<int>().swap(f_mini_pos);
	vector<int>().swap(le);
	vector<int>().swap(re);
	vector<int>().swap(le_r);

	mi = mallinfo2();
	
	double end_ram = mi.hblkhd + mi.uordblks;
	auto end = get_time::now();
	auto diff2 = end - begin;
	output_file << "Construct Time:  "<< chrono::duration_cast<chrono::milliseconds>(diff2).count()<<" ms"<<endl;	
	output_file << "Construct space:" << (end_ram-begin_ram)/1000000 << " MB" << endl;
	
	if(!st.patterns.empty()){
		size_t total_occ_no = 0;
		cout << "pattern matching begin" << endl;
		string pfile = st.patterns;

		ifstream file(pfile, std::ios_base::in | std::ios_base::binary);
		boost::iostreams::filtering_istream patterns;
		patterns.push(boost::iostreams::gzip_decompressor());
		patterns.push(file);	
		begin = get_time::now();
		for (string pattern; getline(patterns, pattern); ){
			if(pattern.size() < k) continue;
			set<int64_t> valid_res;
			// output_file << pattern << ":"; 
			size_t j = find_minimizer_index(pattern, k);
			int l = j;
			int r = pattern.size()-l;
			if(j*2 > pattern.size()){
				string left_pattern = pattern.substr(0, j+1);
				reverse(left_pattern.begin(), left_pattern.end());
				pair<int64_t ,int64_t> left_interval = rev_pattern_matching ( left_pattern, fH, LSA, LLCP, lrmq, (int64_t)g );

				if( left_interval.first <= left_interval.second ){					
					for(int64_t i = left_interval.first; i <= left_interval.second; i++){
						int begin = Nz - (LSA[i]+left_pattern.size());
						if(begin < 0 || begin + pattern.size() >= Nz) continue;
						if(valid_res.count(begin%N)) continue;						
						int c = Nz - LSA[i] - 1;
						double lpi = fH.get_pi(c, begin, l);
						double rpi = fH.check_pi(pattern, l, begin+l, r, c);	
						
						if( lpi * rpi * z >= 1 ){
							valid_res.insert( begin%N );
						}
					}
				}
	
			}else{
				string right_pattern = pattern.substr(j);
				pair<int64_t ,int64_t> right_interval = pattern_matching ( right_pattern, fH, RSA, RLCP, rrmq, (int64_t)g );
				if( right_interval.first <= right_interval.second ){
					for(int64_t i = right_interval.first; i <= right_interval.second; i++){
						int begin = RSA[i] - l;						
						if(begin < 0 || begin + pattern.size() >= Nz) continue;
						if(valid_res.count(begin%N)) continue;
						
						double rpi = fH.get_pi( RSA[i], RSA[i], r);
						double lpi = fH.check_pi(pattern, 0, begin, l,  RSA[i]);
						if( rpi * lpi * z >= 1 ){
							valid_res.insert( begin%N );
						}
					}
				}
			}
			// for(auto occ : valid_res){
				// output_file << occ << " ";
			// }
			// output_file << endl;
			total_occ_no += valid_res.size();
			// return 0;
		}
		end = get_time::now();
		auto diff3 = end - begin;
		output_file << pfile << " Search Time:  "<< chrono::duration_cast<chrono::milliseconds>(diff3).count()<<"ms. \n Totally " << total_occ_no << " occurrences are found." << endl;
	}

	return 0;
}

