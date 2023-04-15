#include <set>
#include <string>
#include <vector>
#include <deque>
#include <iostream>

using namespace std;

int resort(deque<string> const& W, string& m){
	int j = 0;
	m = W[j];
	for(int i = 1; i < W.size(); i++){
		if(W[i] < W[j]){
			m = W[i];
			j = i;
		}
	}
	return j;
}

vector<int> minimizers(string const& S, int k, int w){
	int i = 0;
	deque<string> win_strs;
	set<int> m_pos;
	string minimizer;
	for(int j = 0; j < w; j++){
		string kmer = S.substr(j,k);
		win_strs.emplace_back(kmer);
	}
	int j = resort(win_strs, minimizer);
	m_pos.insert(j);
	for( i = 1; i < S.size()-w-1; i++){
		string kmer = S.substr(i+w-1, k);
		string out = win_strs.front();
		win_strs.pop_front();
		if(out == minimizer){
			j = resort(win_strs,minimizer);
			if(kmer < minimizer){
				m_pos.insert(i+w-1);
				minimizer = kmer;
			}else{
				m_pos.insert(i+j);
			}
		}else{
			if(kmer < minimizer){
				m_pos.insert(i+w-1);
				minimizer = kmer;
			}
		}
		win_strs.emplace_back(kmer);
	}
	return vector<int>(m_pos.begin(), m_pos.end());
}

