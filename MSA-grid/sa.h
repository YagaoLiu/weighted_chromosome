#include <map>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_set>

#include "defs.h"

#include <divsufsort.h>
#include <sdsl/rmq_support.hpp>
#include <boost/functional/hash.hpp>


using namespace std;
using namespace sdsl;

void right_compacted_trie ( unordered_set<int64_t > &anchors, int64_t  * SA, int64_t  * LCP, int64_t  n, int64_t  * RSA, int64_t  * RLCP, int64_t  g );
void left_compacted_trie ( unordered_set<int64_t > &anchors, int64_t  * SA, int64_t  * LCP, int64_t  n, int64_t  * RSA, int64_t  * RLCP, int64_t  g );
int64_t  LCParray ( unsigned char * text, size_t n, int64_t  * SA, int64_t  * ISA, int64_t  * LCP );


