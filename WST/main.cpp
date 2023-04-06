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
#include <chrono>
#include <ctime>

#include "input.h"
#include "weighted_sequence.h"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

using namespace std;
using namespace std::chrono;
using get_time = std::chrono::steady_clock;

int main (int argc, char ** argv ) {
    Settings st = decode_switches(argc, argv);
	istream& text = st.text.is_open()?st.text:cin;
    string& pfile_prefix = st.patterns;
    ostream& output = st.output.is_open()?st.output:cout;
    	
	WeightedSequence W;
    text >> W;

    W.build_index(st.z, st.quiet, output);
	
	cout << "Start patterm matching" << endl;
#if 1
	string pfile_suffix[7] = {"p32.txt.gz","p64.txt.gz","p128.txt.gz","p256.txt.gz","p512.txt.gz","p1024.txt.gz","p2048.txt.gz"};
	for(string ps : pfile_suffix){
		string pfile = pfile_prefix + ps;
		ifstream file(pfile, std::ios_base::in | std::ios_base::binary);
		boost::iostreams::filtering_istream patterns;
		patterns.push(boost::iostreams::gzip_decompressor());
		patterns.push(file);		
		auto begin = get_time::now();
		for (string str; getline(patterns, str); ){
			vector<int> occs = W.occurrences(str);
		}
		auto end = get_time::now();
		auto diff = end - begin;
		output << pfile << " total search:" << chrono::duration_cast<chrono::milliseconds>(diff).count() << endl;
	}
#endif

    return 0;
}

