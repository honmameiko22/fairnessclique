#include<bits/stdc++.h>
using namespace std;
int main(){
	int n, mem;
	vector<vector<int> > arrs;
	ifstream ifs("data/relative_cliques_weak.txt");
	ofstream ofs("data/relative_cliques_filter.txt");
	while(ifs>>n){
		vector<int> arr;
		for(int i=0; i<n; i++) {
			ifs>>mem;
			arr.push_back(mem);
		}
		arrs.push_back(arr);
	}
	for(int i =0 ;i <arrs.size(); i++){
		bool unique=true;
		for(int j=0; j<arrs.size(); j++){
			if(j == i) continue;
			if(includes(arrs[j].begin(),arrs[j].begin()+arrs[j].size(),arrs[i].begin(), arrs[i].begin()+arrs[i].size())){
				if(arrs[i].size()!=arrs[j].size())
					unique = 0;
			}
		}
		if(unique){
			ofs<<arrs[i].size()<<" ";
			for(int j=0; j<arrs[i].size(); j++)
				ofs<<arrs[i][j]<<" ";
			ofs<<"\n";
		}
	}
	return 0;
}