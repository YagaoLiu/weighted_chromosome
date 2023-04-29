/**
 *    Weighted Index
 *    Copyright (C) 2017 Carl Barton, Tomasz Kociumaka, Chang Liu, Solon P. Pissis and Jakub Radoszewski.
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <cstring>
#include <ctime>
#include <chrono>
#include <malloc.h>

#include "input.h"
//#include "weighted_sequence.h"
#include "minimizer_index.h"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

using namespace std;
using namespace std::chrono;

using get_time = std::chrono::steady_clock;

int main (int argc, char ** argv ) {
    Settings st = decode_switches(argc, argv);
    istream& text = st.text.is_open()?st.text:cin;
    ostream& output_file = st.output.is_open()?st.output:cout;
	int ell = st.ell;
	
	auto begin = get_time::now();
	struct mallinfo2 mi;
    mi = mallinfo2();
	auto begin_ram = mi.hblkhd + mi.uordblks;
	MinimizerIndex M;
	text >> M;
    M.build_index(st.z, ell);
	mi = mallinfo2();
	auto end_ram = mi.hblkhd + mi.uordblks;
	auto end = get_time::now();
	auto diff2 = end - begin;
	output_file << "Construct Time:  "<< chrono::duration_cast<chrono::milliseconds>(diff2).count()<<" ms"<<endl;	
	output_file << "Construct space:" << (end_ram-begin_ram)/1000000 << " MB" << endl;

	int total_occ = 0;
	if(!st.patterns.empty()){		
		begin = get_time::now();
		
		ifstream file(st.patterns, std::ios_base::in | std::ios_base::binary);
		boost::iostreams::filtering_istream patterns;
		patterns.push(boost::iostreams::gzip_decompressor());
		patterns.push(file);	
		begin = get_time::now();
		for (string pattern; getline(patterns, pattern); ){
			// output_file << pattern << ": ";
			std::vector<int> occs = M.GridMatch(pattern, ell, st.z);
			// if (occs.empty()) {
				// output_file << "Not found\n";
			// } else {
				// for (auto p : occs) {
					// output_file << p << " ";
				// }
				// output_file << endl;
			// }
			total_occ += occs.size();
		}
		end = get_time::now();
		auto diff = end - begin;
		output_file << "Search Time:  " << chrono::duration_cast<chrono::milliseconds>(diff).count() << "ms. \n Totally " << total_occ << " occurrences are found." << endl;
	}

	

    return 0;
}

