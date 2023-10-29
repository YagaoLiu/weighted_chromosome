#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

int main(){
	string english = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	int length = 10000;
	
	int alphabet_size = 20;
	string alphabet = english.substr(0, alphabet_size);
	ofstream output("random_weighted_string_a20.txt");
	output << length << endl;
	output << alphabet << endl;
	srand(time(NULL));
	
	for(int i = 0; i < length; i++){
		vector<int> prob(alphabet_size, 0);
		int sum = 0;
		int run = 0;
		while(sum < 100){
			run ++;
			int pos = rand()%alphabet_size;
			if(run == 1){
				int q = rand()%21 +80;
				prob[pos] += q;
				sum += q;
			}else if(run == 5){
				prob[pos] += 100 - sum;
				break;
			}else{
				int p = rand()%(101-sum);
				prob[pos] += p;
				sum += p;
			}
		}
		for(int i = 0; i < alphabet_size; i++){
			output << double(prob[i]*0.01) << " ";
		}
		output << endl;
	}
}