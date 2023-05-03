#include <cmath>
#include <chrono>
#include <list>
#include <cstdlib>
#include <ctime>
#include <set>
#include <deque>
#include <unordered_map>
#include <sys/resource.h>
#include <algorithm>
#include "estimation.h"
#include "anchors_new.h"
#include "minimizer_index.h"
#include <unordered_map>

typedef grid_point point;
typedef grid_query query;

using std::endl;
using std::cerr;
using get_time = std::chrono::steady_clock;


bool sort_sa(const pair<int,int> &a,const pair<int,int> &b)
{
       return a.first<b.first;
}

std::istream & operator >> (std::istream& input, MinimizerIndex &M) {
    input >> M.N;
    input >> M.alph;
    int A = M.alph.size();
    for (int i = 0; i < M.N; ++i) {
        double sum = 0;
        std::vector<double> symbol(A, 0);
        for (int j = 0; j < A; ++j) {
            input >> symbol[j];
            sum += symbol[j];
        }
        if (std::abs(sum-1) > EPS) {
            std::cerr << "Probabilities at position " << i << " do not sum up to 1" << std::endl;
            throw 1;
        }
        M.fP.emplace_back(symbol);
    }
 	return input;
}

void extention ( vector<vector<double>>& text, string& s, string& alph, vector<int>& le, vector<int>& re, double z ){
	int n = text.size();
	int N = s.size();
	unordered_map<char, int> A;
	for(int i = 0; i < alph.size(); i++){
		A[alph[i]] = i;
	}
	
	vector<double> _pi;
	for(int i = 0; i < N; i++){
		_pi.push_back( text[i%n][A[s[i]]]);
	}	

	re.assign(N,0);	
	for (auto j = 0; j < z; j++ ){
		int i = 0;
		int si = i + j*n;
		int l = 0;
		double cum_pi = _pi[si];
		while( (cum_pi * z >= 1) && (i+l < n) ){
			l++;
			cum_pi *= _pi[si+l];
		}
		re[si] = l-1<0 ? 0 : l-1;
		
		for(i = 1; i < n; i++){
			si++;
			if(l > 0) l--;
			cum_pi /= _pi[si-1];
			while( (cum_pi * z >= 1) && (i+l < n) ){
				l++;
				cum_pi *= _pi[si+l];
			}
			re[si] = l-1<0 ? 0 : l-1;
		}		
	}
	
	le.assign(N,0);
	reverse(_pi.begin(), _pi.end());
	for (auto j = 0; j < z; j++ ){
		int i = 0;
		int si = i + j*n;
		int l = 0;
		double cum_pi = _pi[si];
		while( (cum_pi * z >= 1) && (i+l < n) ){
			l++;
			cum_pi *= _pi[si+l];
		}
		le[si] = l-1<0 ? 0 : l-1;
		
		for(i = 1; i < n; i++){
			si++;
			if(l > 0) l--;
			cum_pi /= _pi[si-1];
			while( (cum_pi * z >= 1) && (i+l < n) ){
				l++;
				cum_pi *= _pi[si+l];
			}
			le[si] = l-1<0 ? 0 : l-1;
		}		
	}
	reverse(le.begin(), le.end());
}

void MinimizerIndex::build_index(double z, int ell){
	vector<vector<double>> rP(fP.rbegin(), fP.rend());
	Estimation fS(fP,alph,z);
	PropertyString fT;
	std::vector<int> f_mini_pos;
	std::vector<int> r_mini_pos;
	int i = 0;
	int k = ceil(4 * log2(ell) / log2(alph.size()));
	int w = ell - k + 1;
	
	for(PropertyString const & s : fS.strings()){
		fT += s;
		std::vector<int> M;
		minimizers_with_kr(s.string(), M, w, k);		
		for(auto it : M){
			f_mini_pos.emplace_back(it + i*N);
		}
		i++;
	}	
	string zstrs = fT.string();
	string rev_zstrs(zstrs.rbegin(), zstrs.rend());
	vector<int> le;
	vector<int> re;
	extention(fP, zstrs, alph, le, re, z);
	vector<int> le_r(le.rbegin(), le.rend());
	vector<int> re_r(re.rbegin(), re.rend());

	f_mini_pos.erase(remove_if(f_mini_pos.begin(), f_mini_pos.end(), [&](int i) { return re[i] < k; }), f_mini_pos.end());

	for(auto i : f_mini_pos){
		r_mini_pos.push_back(zstrs.size() - i - 1);
	}
	
	HeavyString fH(fP, zstrs, alph, f_mini_pos, le, re, true);
	HeavyString rH(rP, rev_zstrs, alph, r_mini_pos, le_r, re_r, false);
	Nz = zstrs.size();
	int   g = f_mini_pos.size();
	
	forward_index = new PropertySuffixTree(re, fH, f_mini_pos);	
	reverse_index = new PropertySuffixTree(le_r, rH, r_mini_pos);
	
	RSA = forward_index->toSA();
	vector<int> LSA = reverse_index->toSA();
		
	vector<pair<int  ,int  >> l;
	vector<pair<int  ,int  >> r;
	for ( int   i = 0; i < g; i++ )
  	{
  		l.push_back( make_pair( Nz-1-LSA[i], i) );
 		r.push_back( make_pair( RSA[i], i ) );	
	}
	
 	sort(l.begin(),l.end(),sort_sa);
	sort(r.begin(),r.end(),sort_sa);

  	std::vector<point> points;
  	  	
  	for ( int   i = 0; i < g; i++ )
  	{
 
		point to_insert;
		to_insert.row = l.at(i).second+1;
		to_insert.col = r.at(i).second+1;
		
		to_insert.level = 0;
		to_insert.label = l.at(i).first;
		points.push_back(to_insert); 

  	}
  	  	
  	//constructing grid with bd-anchors
  	construct.build( points, 0 );

	fS.clear();
	fT.clear();	
	vector<vector<double>>().swap(fP);
	
		int oc = zstrs.find("AGCTAGCTTCAAATAGTATTCTGCAGCCAGGC");
	cout << oc << " " << oc%N <<endl;
	cout << fH.substr(oc, 32) << endl;
}

int  find_minimizer_index(string s, int  k) {
  int  minimizer_index = 0;
  for (int  i = 1; i <= s.length() - k; i++) {
    if (s.substr(i, k) < s.substr(minimizer_index, k)) {
      minimizer_index = i;
    }
  }

  return minimizer_index;
}

std::vector<int> MinimizerIndex::occurrences(std::string const &P, int ell, double z, std::ostream& result) const{
	int m = P.size();
	int k = ceil(4 * log2(ell) / log2(alph.size()));
	std::string minimizer = P.substr(0,k);
	int min_index = 0;
	for(size_t i = 1; i < m-k; i++){
		std::string kmer = P.substr(i,k);
		if(kmer < minimizer){
			minimizer = kmer;
			min_index = i;
		}
	}
	std::vector<int> final_occs;
	std::set<int> occs;
	if(min_index <= ceil(m/2) ){
		std::string min_suf = P.substr(min_index);
		for(int o : forward_index->occurrences(min_suf)){
			if( (o%N) >= min_index ){
				occs.insert( (o%N) - min_index );
			}
		}
	}else{
		std::string rP = P;
		std::reverse(rP.begin(), rP.end());
		int r_index = m - min_index - 1;
		std::string min_suf = rP.substr(r_index);
		for(int o : reverse_index->occurrences(min_suf)){
			if( (o%N) >= r_index ){
				occs.insert( (o%N) - r_index );
			}
		}
	}
	for(int i : occs){
		if(is_valid(P, i, z)){
			final_occs.emplace_back(i);
		}
	}
	return final_occs;
}


std::vector<int> MinimizerIndex::GridMatch(std::string const& pattern, int ell, double z) const{
	int m = pattern.size();
	int k = ceil(4 * log2(ell) / log2(alph.size()));
	std::string minimizer = pattern.substr(0,k);
	int min_index = find_minimizer_index(pattern,k);
	std::set<int> occs;
	std::string left_pattern = pattern.substr(0, min_index+1);
	reverse(left_pattern.begin(), left_pattern.end());
	pair<int,int> left_interval = reverse_index->SAoccurrences(left_pattern);
	string right_pattern = pattern.substr(min_index);
	pair<int,int> right_interval = forward_index->SAoccurrences(right_pattern);

	if ( left_interval.first <= left_interval.second  && right_interval.first <= right_interval.second )
	{	
		// Finding rectangle containing bd-anchors in grid
		grid_query rectangle;
		rectangle.row1 = left_interval.first+1;
		rectangle.row2 = left_interval.second+1;
		rectangle.col1 = right_interval.first+1;
		rectangle.col2 = right_interval.second+1;			
		
		vector<size_t> result;
		construct.search_2d(rectangle, result);
		
		set<int64_t> valid_res;	
		for(size_t i = 0; i < result.size(); i++){
			if(RSA[result.at(i)-1] < min_index) continue;
			if(occs.count((RSA[result.at(i)-1]-min_index)%N)) continue;
			double pi = forward_index->get_pi(RSA[result.at(i)-1], RSA[result.at(i)-1]-min_index, pattern.size());
			if(pi * z >= 1){
				occs.insert((RSA[result.at(i)-1]-min_index)%N);
			}
		}
	}
	std::vector<int> final_occs(occs.begin(), occs.end());
	return final_occs;
}

bool  MinimizerIndex::is_valid(string const& p, int pos, double z) const{
		unordered_map<char, int> mapping;
		pos = pos%N;
		if(pos + p.size() > N) return false;
		for(int i = 0; i < alph.size(); i++){
			mapping[alph[i]] = i;
		}
		double cum_prob = 1;
		for(int i = 0; i < p.size(); i++){
			cum_prob *= fP[pos+i][mapping[p[i]]];
		}
		return (cum_prob*z >= 1)?true:false;
}