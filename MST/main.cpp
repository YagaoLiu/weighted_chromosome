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

	begin = get_time::now();
	if(st.patterns.is_open()){
		istream& pattern_file = st.patterns;
		string pattern;
		while(true){
			if (!(pattern_file >> pattern)) break;
			//output << pattern << ": ";
			std::vector<int> occs = M.occurrences(pattern, ell, st.z, output_file);
			// if (occs.empty()) {
				// output << "Not found\n";
			// } else {
				// for (auto p : occs) {
					// output << p << " ";
				// }
				// output << endl;
			// }
		}
	}
	end = get_time::now();
	auto diff = end - begin;
	output_file << "Search Time:  " << chrono::duration_cast<chrono::milliseconds>(diff).count() << endl;
	

    return 0;
}

