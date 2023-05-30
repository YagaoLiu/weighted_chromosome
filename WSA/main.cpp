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
#include <malloc.h>

#include "node.h"
#include "defs.h"
#include "input.h"
#include "estimation.h"
#include "property_string.h"
#include "match.h"

#include <divsufsort.h>
#include <sdsl/rmq_support.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

using namespace std;
using namespace std::chrono;
using get_time = chrono::steady_clock;

int main (int argc, char ** argv )
{
    Settings st = decode_switches(argc, argv);
	istream& text_file = st.text.is_open()?st.text:cin;
	ostream& output_file = st.output.is_open()?st.output:cout;
    string& pfile_prefix = st.patterns;
	ofstream result;
	
	double z = st.z;
	
	auto begin = get_time::now();
	struct mallinfo2 mi;
    mi = mallinfo2();
	double begin_ram = mi.hblkhd + mi.uordblks;
	string alphabet;
	vector<vector<double>> text;
	
	size_t N;
	text_file >> N;
	text_file >> alphabet;
	size_t A = alphabet.size();
	for (size_t i = 0; i < N; ++i) {
        double sum = 0;
        vector<double> symbol(A, 0);
        for (size_t j = 0; j < A; ++j) {
            text_file >> symbol[j];
            sum += symbol[j];
        }
        text.emplace_back(symbol);
    }

	// cout << "text length:" << N << endl;
	// cout << "Finish reading input" << endl;
	

	// cout << "index begin" << endl;
	Estimation fS(text, alphabet, z);
	PropertyString fT;
	
	for(PropertyString const & s : fS.strings()){
		fT += s;
	}
	
	fS.clear();
	string zstrs = fT.string();
	
	size_t Nz = fT.string().size();

	int * SA	= new int [Nz];
	int * LCP	= new int [Nz];
	int * ME	= fT.property_pointer();
		
	begin = get_time::now();
	SA_LCP_index ( text, fT.string().c_str(), Nz, N, z, SA, LCP );
//	union_find_resort ( SA, LCP, ME, Nz );

	vector<int> tmp_lcp ( LCP, LCP+Nz);
	rmq_succinct_sct<> rmq ( &tmp_lcp );
	vector<vector<double>>().swap(text);
	
	mi = mallinfo2();
	
	double end_ram = mi.hblkhd + mi.uordblks;
	auto end = get_time::now();
	auto diff2 = end - begin;
	output_file<<"CT "<< chrono::duration_cast<chrono::milliseconds>(diff2).count()<< endl;	
	output_file << "IS " << (end_ram-begin_ram)/1000000 << endl;
	
	string pfile_suffix[7] = {"p32.txt.gz","p64.txt.gz","p128.txt.gz","p256.txt.gz","p512.txt.gz","p1024.txt.gz","p2048.txt.gz"};
	for(string ps : pfile_suffix){	
		string pfile = pfile_prefix + ps;
		begin = get_time::now();
		ifstream file(pfile, std::ios_base::in | std::ios_base::binary);
		boost::iostreams::filtering_istream patterns;
		patterns.push(boost::iostreams::gzip_decompressor());
		patterns.push(file);		
		auto begin1 = get_time::now();
		for (string pattern; getline(patterns, pattern); ){
			int m = pattern.size();
			list<int> Occ;
			//output_file << pattern << ":\n";
			int num_Occ = match ( pattern, zstrs, N, SA, LCP, ME, Occ, rmq );
			if ( num_Occ == 0 )	{
				//output_file << "No found\n";
			}else{
				for ( auto it = Occ.begin(); it != Occ.end(); it++ )
				{
					if(ME[*it] >= pattern.size()){
						//output_file << *it << ' ';
					}
				}
				//output_file << '\n';
			}
			Occ.clear();
		}
		end = get_time::now();
		auto diff3 = end - begin;
		output_file << "PMT "<< chrono::duration_cast<chrono::milliseconds>(diff3).count()<< endl;
	}	

#if 0 
	vector<int> WSA;
	vector<int> WLCP;
	WSA.reserve ( N );
	WLCP.reserve ( N );
	begin =get_time:: now();
	WeightedIndex ( n, N, SA, LCP, ME, WSA, WLCP );
	

	for ( int i = 0; i < N; i++ )
		if ( ME[SA[i]] != 0 )
			result << fT.string().substr(SA[i], ME[SA[i]]) << '\n';
	result.close();	
#endif
	
	return 0;
}

