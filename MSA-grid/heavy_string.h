#ifndef HEAVY_STRING_H 
#define HEAVY_STRING_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <stdexcept> 
#include <iostream>

class HeavyString{
	std::string H;
	std::unordered_map<int, char> _alt;
	std::unordered_map<int, double> delta_pi;
	std::unordered_map<int, std::vector<int>> alt_pos;
	std::vector<double> pi_suf;
	int n;
	int N;

	public:
	HeavyString(): n(0), N(0){};
	
	HeavyString(std::vector<std::vector<double>>& P, std::string const& S, std::string& A) 
		: n(P.size()), N(S.size()) {
			if (n == 0 || N == 0) {
				throw std::invalid_argument("P and S cannot be empty.");
			}

			for(size_t i = 0; i < n; i++){
				int which_max = max_element(P[i].begin(), P[i].end()) - P[i].begin();
				H+=(A[which_max]);
			}		
			for(size_t i = 0; i < N; i++){
				if(S[i] != H[i%n]){
					if (A.find(S[i]) == std::string::npos) {
						throw std::invalid_argument("S contains a character not in A.");
					}
					_alt[i] = S[i];
				}
			}
		}	
		
	HeavyString(std::vector<std::vector<double>>& P,  std::string const& S, std::string& A, std::unordered_set<int> min_pos, std::vector<int> le, std::vector<int> re){
		n = P.size();
		N = S.size();
		if (n == 0 || N == 0) {
			throw std::invalid_argument("P and S cannot be empty.");
		}
		
		std::vector<double> pi_arr;
		
		for(size_t i = 0; i < n; i++){
			int which_max = max_element(P[i].begin(), P[i].end()) - P[i].begin();
			H+=(A[which_max]);
			double pi = log2(P[i][which_max]);
			pi_arr.push_back(pi);
		}
		
		
		pi_suf.assign(pi_arr.begin(), pi_arr.end());
		for(int i = n-2; i >= 0; i--){
			pi_suf[i] += pi_suf[i+1];
		}
		for(size_t i = 0; i < N; i++){
			if(S[i] != H[i%n]){
				if (A.find(S[i]) == std::string::npos) {
					throw std::invalid_argument("S contains a character not in A.");
				}
				_alt[i] = S[i];
			}
		}
		for(int m : min_pos){
			int begin = m - le[m];
			int end = m + re[m] + 1;
		//	std::cout << m << " " << begin << " " << end << std::endl;
			for(int i = begin; i < end; i++){
				int h = i%n;
				if(H[h] != S[i]){
					double this_pi = log2(P[h][A.find(S[i])]);
					//_alt[i] = S[i];
					alt_pos[m].push_back(i);
					delta_pi[i] =  this_pi - pi_arr[i];
				}
			}
		}
	}
	
	HeavyString(const HeavyString& other): H(other.H), _alt(other._alt), n(other.n), N(other.N) {}

	HeavyString& operator=(const HeavyString& other) {
		if (this != &other) {
			H = other.H;
			n = other.n;
			N = other.N;
			std::unordered_map<int, char> temp(other._alt);
			std::swap(_alt, temp);
		}
		return *this;
	}

	char& operator[](size_t i){ 
		if (i >= N) {
			throw std::out_of_range("Index out of range.");
		}

		if(_alt.count(i)){
			return _alt.at(i);
		}else{
			return H[i%n];
		}
	}
	
	std::string substr(size_t pos){
		if(pos >= N){
			return "";
		}
		std::string substring = H.substr(pos%n);
		for(size_t i = 0; i < substring.size(); i++){
			if(_alt.count(pos+i)){
				substring[i] = _alt.at(pos+i);
			}
		}
		return substring;			
	}

	std::string substr(size_t pos, size_t len = 0){
		if (pos >= N || pos%n + len > n) {
			return "";
		}
		if(len == 0){	//if no length is given, print to the end
			len = n- pos%n;
		}
		std::string substring = H.substr(pos%n, len);

		for(size_t i = 0; i < len; i++){
			if(_alt.count(pos+i)){
				substring[i] = _alt.at(pos+i);
			}
		}

		return substring;
	}
	
	double get_pi(int i, int begin, int length){
		if(begin%n > i%n)
			return 0;
		if(begin%n + length > n)
			return 0;
		
		
		int end = begin%n + length;
		
		double cum_pi = pi_suf[begin%n] - pi_suf[begin%n + length];
		
		std::vector<int>& v = alt_pos[i];
		if(v.empty()){
			return pow(2,cum_pi);
		}else{
			auto lower = std::lower_bound(v.begin(), v.end(), begin);
			auto upper = std::upper_bound(v.begin(), v.end(), end);		
			
			if(lower < upper){
				for(auto it = lower; it != upper; it++){
					cum_pi += delta_pi[*it];
				}
			}
			return pow(2,cum_pi);
		}
	}

	size_t length() const {return N;}
	size_t size() const {return N;}
	size_t heavy_length() const {return n;}

	class Iterator : public std::iterator<std::random_access_iterator_tag, char> {
		HeavyString* hs; 
		size_t index;
		public:
		Iterator(HeavyString* hs, size_t index) : hs(hs), index(index) {}

		char& operator*() {
			return (*hs)[index];
		}

		Iterator& operator++() { 
			++index;
			return *this;
		}

		Iterator operator++(int) {
			Iterator temp = *this;
			++index;
			return temp;
		}

		Iterator& operator--() {
			--index;
			return *this;
		}

		Iterator operator--(int) {
			Iterator temp = *this;
			--index;
			return temp;
		}

		Iterator operator+(size_t n) const {
			return Iterator(hs, index + n);
		}

		friend Iterator operator+(size_t n, const Iterator& it) {
			return it + n;
		}

		Iterator& operator+=(size_t n) {
			index += n;
			return *this;
		}

		Iterator operator-(size_t n) const {
			return (*this) + (-n);
		}

		friend ptrdiff_t operator-(const Iterator& lhs, const Iterator& rhs) {
			return lhs.index - rhs.index; 
		}

		Iterator& operator-=(size_t n) {
			(*this) += (-n);
			return *this;
		}

		bool operator<(const Iterator& other) const {
			return index < other.index; 
		}

		bool operator>(const Iterator& other) const {
			return other < (*this); 
		}

		bool operator<=(const Iterator& other) const {
			return !((*this) > other); 
		}

		bool operator>=(const Iterator& other) const {
			return !((*this) < other); 
		}

		bool operator==(const Iterator& other) const {
			return index == other.index; 
		}

		bool operator!=(const Iterator& other) const {
			return !((*this)==other); 
		}   
	};


	Iterator begin() {return Iterator(this, 0);}
	Iterator end() {return Iterator(this,N);}
};


std::ostream& operator<<(std::ostream& os, HeavyString& hs);

#endif