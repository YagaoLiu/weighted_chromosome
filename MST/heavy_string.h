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
	std::unordered_map<size_t, char> _alt;
	size_t n;
	size_t N;

	public:
	HeavyString(std::vector<std::vector<double>>& P, std::string const& S, std::string& A) 
		: n(P.size()), N(S.size()) {
			if (n == 0 || N == 0) {
				throw std::invalid_argument("P and S cannot be empty.");
			}

			for(size_t i = 0; i < n; i++){
				size_t which_max = max_element(P[i].begin(), P[i].end()) - P[i].begin();
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
	
	HeavyString(std::string const& S): n(S.size()), N(S.size()){
			H = S;
	}

	HeavyString(const HeavyString& other): H(other.H), _alt(other._alt), n(other.n), N(other.N) {}

	HeavyString& operator=(const HeavyString& other) {
		if (this != &other) {
			H = other.H;
			n = other.n;
			N = other.N;
			std::unordered_map<size_t, char> temp(other._alt);
			std::swap(_alt, temp);
		}
		return *this;
	}

	char& operator[](size_t i) { 
		if (i >= N) {
			throw std::out_of_range("Index out of range.");
		}

		if(_alt.count(i)){
			return _alt.at(i);
		}else{
			return H[i%n];
		}
	}

	std::string substr(size_t pos, size_t len){
		if (pos >= N || len == 0 || pos + len > n) {
			return "";
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
	size_t heavy_length() const {return n;}

	class const_iterator : public std::iterator<std::random_access_iterator_tag, char> {
    private:
        size_t index;
        HeavyString* hs;
    public:
        // constructor
		const_iterator(size_t i=-1, HeavyString* h=nullptr) : index(i), hs(h) {}
		
        // dereference operator
        char operator*() const  {
            return (*hs)[index];
        }

        // increment operators
        const_iterator& operator++() {
            index++;
            return *this;
        }

        const_iterator operator++(int) {
            const_iterator temp = *this;
            ++(*this);
            return temp;
        }

        // decrement operators
        const_iterator& operator--() {
            index--;
            return *this;
        }

        const_iterator operator--(int) {
            const_iterator temp = *this;
            --(*this);
            return temp;
        }

		const_iterator operator+(size_t n) const {
			return const_iterator(index + n, hs);
		}

		friend const_iterator operator+(size_t n, const const_iterator& it) {
			return it + n;
		}

		const_iterator& operator+=(size_t n) {
			index += n;
			return *this;
		}

		const_iterator operator-(size_t n) const {
			return (*this) + (-n);
		}

		friend ptrdiff_t operator-(const const_iterator& lhs, const const_iterator& rhs) {
			return lhs.index - rhs.index; 
		}

		const_iterator& operator-=(size_t n) {
			(*this) += (-n);
			return *this;
		}

        // comparison operators
        bool operator==(const const_iterator& other) const  {
            return (index == other.index && hs == other.hs);
        }

       bool operator==(const std::string::const_iterator& other) const  {
            return ((*hs)[index] == *other);
        }
		
        bool operator!=(const const_iterator& other) {
            return !(*this == other);
        }
	};

	const_iterator begin() {return const_iterator(0, this);}
	const_iterator end() {return const_iterator(N, this);}
};

std::ostream& operator<<(std::ostream& os, HeavyString& hs);
#endif