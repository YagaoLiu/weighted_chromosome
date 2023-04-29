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

void property_computer(vector<vector<double>>& matrix, string& text, string& A, vector<int>& _pi, double z){
	int n = matrix.size();
	int R = text.size() / n;
	unordered_map<char, int> a;
	for(size_t i = 0; i < A.size(); i++){
		a[A[i]] = i;
	}
	
	for(int r = 0; r < R; r++){
		double product = 1.0;
		deque<double> q;
		vector<int> _p(n,0);
		string subtext = text.substr(r*n,n);
		for(int i = 0; i < n; i++){
			if(i == 0){
				int j = i;				
				while(  (j < n) && (product * matrix[j][a[subtext[j]]] * z >= 1) ){
					product *= matrix[j][a[subtext[j]]];
					q.push_back(matrix[j][a[subtext[j]]]);
					j++;
				}
				_p[i] = q.size();
			}else if(i+q.size() > n){
				_p[i] = _p[i-1] - 1;
			}else{
				if(!q.empty()){
					double last = q.front();
					product /= last;
					q.pop_front();
				}
				int j = i + q.size();
				while(  (j < n) && (product * matrix[j][a[subtext[j]]] * z >= 1) ){
					product *= matrix[j][a[subtext[j]]];
					q.push_back(matrix[j][a[subtext[j]]]);
					j++;
				}
				_p[i] = q.size();			
			}
		}
		_pi.insert(_pi.end(), _p.begin(), _p.end());
	}					
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
		re[si] = l-1;
		
		for(i = 1; i < n; i++){
			si++;
			l--;
			cum_pi /= _pi[si-1];
			while( (cum_pi * z >= 1) && (i+l < n) ){
				l++;
				cum_pi *= _pi[si+l];
			}
			re[si] = l-1;
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
	for(auto i : f_mini_pos){
		r_mini_pos.push_back(zstrs.size() - i);
	}
	// property_computer(fP, zstrs, alph, f_pi,z);
	// property_computer(rP, rev_zstrs, alph, r_pi,z);
	extention(fP, zstrs, alph, le, re, z);
	vector<int> le_r(le.rbegin(), le.rend());
	vector<int> re_r(re.rbegin(), re.rend());
	
	HeavyString fH(fP, zstrs, alph, f_mini_pos, le, re, true);
	HeavyString rH(rP, rev_zstrs, alph, r_mini_pos, le_r, re_r, false);
	Nz = zstrs.size();
	
	forward_index = new PropertySuffixTree(re, fH,f_mini_pos);
	
	reverse_index = new PropertySuffixTree(le_r, rH, r_mini_pos);
		
	fS.clear();
	fT.clear();	
	vector<vector<double>>().swap(fP);
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
	int min_index = find_minimizer_index(P,k);
	
	
	std::set<int> occs;
	int l = min_index;
	int r = m - min_index;
	if(min_index <= ceil(m/2) ){
		std::string min_suf = P.substr(min_index);

		for(int o : forward_index->occurrences(min_suf)){
			int begin = o - l;
			if(occs.count(begin%N)) continue;
			double rpi = forward_index->get_pi(o, o, r);
			double lpi = forward_index->naive_check(P, 0, begin, l, o);
			
			if( rpi * lpi *z >=1 ){
				occs.insert( begin%N );
			}
		}
	}else{
		std::string rP = P;
		std::reverse(rP.begin(), rP.end());
		int r_index = m - min_index;
		std::string min_suf = rP.substr(r_index);
		for(int o : reverse_index->occurrences(min_suf)){
			int begin = Nz - ( o + min_suf.size() - 1 ) - 1;
			if(occs.count(begin%N)) continue;
			int c = Nz - o;
			double lpi = forward_index->get_pi(c,begin,l);
			double rpi = forward_index->naive_check(P, l, begin+l, r, c);			
			if( rpi * lpi *z >=1 ){
				occs.insert( begin%N );
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