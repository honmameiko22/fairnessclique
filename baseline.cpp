#include<iostream>
#include<set>
#include "LinearHeap.h"

#include "Timer.h"
#include "Graph.h"
#include "Utility.h"
using namespace std;

set<int> intersection(set<int> set_a, set<int> set_b) {
    set<int> result;
    set_intersection(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(), inserter(result, result.begin()));
    return result;
}

void printset(set<int> s){
	for(auto it = s.begin(); it!=s.end(); it++){
		cout<<*it<<" ";
	}
	cout<<endl;
}

vector<int> Graph::degeneracy_ordering(vector<int> vertices){
	printf("compute degeneracy...\n");
	int* ordered_c = new int[n];
	int* degree_arr = new int[n];
	int* vis = new int[n];

	for(int i=0; i<n; i++){
		ordered_c[i] = i;
		degree_arr[i] = 0;
		vis[i]= 0;
	}
	for(int i=0; i<vertices.size(); i++){
		degree_arr[vertices[i]] = pend[vertices[i]] - offset[vertices[i]];
		vis[vertices[i]]=1;
	}

	ListLinearHeap *heap = new ListLinearHeap(n, n);
	heap->init(n, n, ordered_c, degree_arr);
	vector<int> peeling_order;
	peeling_idx = new int[n];
	int cnt=0;
	for(int i=0; i<n; i++){
			int nod, key;
			heap->pop_min(nod, key);
			if(vis[nod]){
				vis[nod]=0;
				peeling_order.push_back(nod);
				peeling_idx[nod]= cnt++;
				//printf("peeling idx[%d] =%d\n", nod, peeling_idx[nod]);
				for(int j=offset[nod]; j<pend[nod]; j++){
						if(vis[edge_list[j]]==1){
							heap->decrement(edge_list[j], 1);
						}
				}
			}
	}
	//printf("degeneracy ordering = ");
	//for(int i=0; i<peeling_order.size(); i++)
	//	printf("%d ", peeling_order[i]);
	//printf("\n");
	delete[] ordered_c;
	delete[] degree_arr;
	delete[] vis;
	return peeling_order;
}

void Graph::BronKerbosch(set<int>R, set<int>P, set<int>X, vector<set<int> >& mc){
	change_mem(sizeof(int)*R.size(),1);
	change_mem(sizeof(int)*P.size(),1);
	change_mem(sizeof(int)*X.size(),1);
	/*
	printf("R:");
	printset(R);
	printf("P:");
	printset(P);
	*/
	if(P.size()==0 && X.size()==0){
		if(R.size()>0){
			int* cnt = new int[attr_size];
			for(int i=0; i<attr_size; i++) cnt[i] = 0;
			for(auto it=R.begin(); it!=R.end(); it++){
				cnt[attribute[*it]]++;
			}
			for(int i=0; i<attr_size; i++){
				if(cnt[i]<threshold) return;
			}
			mc.push_back(R);
		}
		return ;
	}
	
	int maxnode=-1, maxcnt=-1;
	for(auto node : P){
		int localcnt=0;
		for(int i=offset[node]; i<pend[node]; i++){
			if(P.count(edge_list[i])){
				localcnt++;
			}
		}
		if(localcnt > maxcnt){
			maxcnt = localcnt; maxnode = node;
		}
	}
	if(maxnode < 0) return;
	set<int> not_visit;
	
	for(int i= offset[maxnode]; i<pend[maxnode]; i++){
		not_visit.insert(edge_list[i]);
	}
	change_mem(sizeof(int)*not_visit.size(),1);
	
	
	set<int> P_dup =P;
	for(auto it = P_dup.begin(); it!=P_dup.end(); it++){
		int me =*it;
		if(not_visit.count(me)) continue;
		if(!P.count(me))
			continue;
		set<int> my_nei;
    for(int j=offset[me]; j<pend[me]; j++){
      my_nei.insert(edge_list[j]);
    }
		change_mem(sizeof(int)*my_nei.size(),1);
		R.insert(me);
		BronKerbosch(R, intersection(P_dup, my_nei), intersection(X, my_nei), mc);
		R.erase(me);
		P.erase(me);
		X.insert(me);
		change_mem(sizeof(int)*my_nei.size(),0);
	}
	change_mem(sizeof(int)*not_visit.size(),0);
	change_mem(sizeof(int)*R.size(),0);
	change_mem(sizeof(int)*P.size(),0);
	change_mem(sizeof(int)*X.size(),0);
}


bool Graph::checksingle(set<int> clique){
	int* count = new int[attr_size+1];
	for(int i =0 ;i<attr_size+1 ; i++)
		count[i] = 0;
	for(auto c:clique){
		count[attribute[c]]++;
	}
	for(int i=0; i<attr_size; i++){
		if(count[i]<threshold){
			delete[] count;
			return false;
		}
	}
	delete[] count;
	return true;
}
void Graph::Baseline(const char* algo, int delta ){
	vector< set<int> > maximal_cliques;
	printf("threshold = %d, delta=%d\n", threshold,delta);
	Timer t;
	//shrink graph
	colorfulcore();
	//get degeneracy ordering
	//degeneracy_ordering(left);
	//bk-algorithm
	set<int> R, P, X;
	vector<int> vR, vX;
	for(int i=0; i<left.size(); i++){
		P.insert(left[i]);
	}
	BronKerbosch(R,P,X,maximal_cliques);
	printf("number of fair clique =%d\n", maximal_cliques.size());
	/*
	string tmp_path = mid_path;
	ofstream hfs(tmp_path);
	for(int i=0; i<maximal_cliques.size(); i++){
		for(auto c:maximal_cliques[i]){
			hfs<<c<<" ";
		}
		hfs<<endl;
	}
*/
	Timer t1;
	//filter
	if(strcmp(algo,"weak")==0){
		for(int i=0; i<maximal_cliques.size(); i++){
			if(checksingle(maximal_cliques[i])){
				vector<int> vC (maximal_cliques[i].begin(), maximal_cliques[i].end());
				ResultSet.push_back(vC);
			}
		}
	}else if(strcmp(algo,"strong")==0){
		printf("strong\n");
 		nvis=new int[n];
    rvis=new int[n];
    X_vis=new int[n];
    int* vis = new int[n];
    nei_vis=new int[n];
    for(int i=0; i<n; i++) X_vis[i]=0;
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<n; i++) nvis[i]=0;
    for(int i=0; i<n; i++) rvis[i]=0;
    for(int i=0; i<n; i++) nei_vis[i]=0;
    for(int i=0; i<left.size(); i++)
        vis[left[i]]=1;
    ResultSet.clear();

		for(int i=0; i<maximal_cliques.size(); i++){
			vector<int> candidates[attr_size];
			for(auto it = maximal_cliques[i].begin(); it !=maximal_cliques[i].end(); it++){
					candidates[attribute[*it]].push_back(*it);
			}
			Backtrack_Strong(vR, candidates, vX, 0);
		}
		printf("tot time for process: %s (microseconds)\n", Utility::integer_to_string(t1.elapsed()).c_str());

		delete[] rvis;
    delete[] nvis;
    delete[] vis;
    delete[] X_vis;
	}else if(strcmp(algo,"relative")==0){
		ofstream relative_ofs("data/relative_cliques_weak.txt");
		printf("max mem =%d\n", max_mem);
		printf("start finding relative-weak fair cliques\n");
		printf("size of weak fair cliques is %d\n", maximal_cliques.size());
		int relative_num = 0;
    assert(maximal_cliques.size()!=0);
    for(int i=0; i<maximal_cliques.size(); i++){
				vector<int> clique;
				for(set<int>::iterator it = maximal_cliques[i].begin(); it!=maximal_cliques[i].end(); it++)
					clique.push_back(*it);
        enum_related_w(clique, delta, relative_ofs);
    }
		printf("relative weak num =%d\n", relative_num);
		printf("time for relative fair cliques(using weak): %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
		relative_ofs.close();
	}

	printf("tot time for search: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
}
