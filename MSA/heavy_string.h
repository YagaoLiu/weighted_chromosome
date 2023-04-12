#ifndef HEAVY_STRING_H 
#define HEAVY_STRING_H

#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstddef>
#include <stdexcept> 
#include <iostream>

class HeavyString{
	std::string H;
	std::unordered_map<int, char> _alt;
	int n;
	int N;

	public:
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