#include "Utility.h"
#include "Graph.h"
#include "LinearHeap.h"
#include "Timer.h"
using namespace std;

int relative_num;

void Graph::VerifyRClique(const char* attrname, const char* cliques,int attr_size, int thresh){
	printf("start\n");
	ifstream attr_ifs(attrname);
	ifstream clique_ifs(cliques);
	int a, b;
	map<int, int >attr_map;
	while(attr_ifs>>a>>b){
		attr_map[a]=b;
	}
	int num, node;
	int* attr_cnt = new int[attr_size];
	set<set<int> > compare;
	while(clique_ifs>>num){
		vector<int> cur_clique;
		set<int> sorted_clique;
		for(int i=0; i<num; i++) attr_cnt[i]=0;
		for(int i=0; i<num; i++){
			clique_ifs>>node;
			cur_clique.push_back(node);
			sorted_clique.insert(node);
			attr_cnt[attr_map[node]]++;
		}
		if(compare.find(sorted_clique)!=compare.end()) printf("dup!\n");
		compare.insert(sorted_clique);
		int min_=100, max_=-1;
		for(int i=0; i<attr_size; i++){
			if(attr_cnt[i]<min_) min_ = attr_cnt[i];
			if(attr_cnt[i]>max_) max_ = attr_cnt[i];
		}
		if(max_- min_ > thresh) {
			printf("error! sub =%d\n", max_-min_);
			printf("clique: ");
			for(auto i:cur_clique){
				printf("%d(%d) ",i, attr_map[i]);
			}
			printf("\n");
		}
	}
	printf("verify done\n");
}

void Graph::deep_search(vector<int> clique, vector<int>& enum_attrs, vector<vector<int>>& enum_nodes, int ptr, int pos, int limit, ofstream& relative_ofs, int depth){
	if(ptr == enum_attrs.size()) {
        if(clique.size()<threshold*attr_size) return;
        relative_num ++;
        relative_ofs<<clique.size()<<" ";
        printClique(clique, relative_ofs);
        //for(int i=0; i<clique.size(); i++)
			//		relative_ofs<<clique[i]<<" ";
			//	relative_ofs<<"\n";
        return;
    }
    for(int i=pos; i<enum_nodes[ptr].size(); i++){
        clique.push_back(enum_nodes[ptr][i]);
        if(depth+1 < limit)
            deep_search(clique, enum_attrs, enum_nodes, ptr, i+1, limit,relative_ofs, depth+1);
        else deep_search(clique, enum_attrs, enum_nodes, ptr+1, 0, limit, relative_ofs, 0);
		clique.pop_back();
    }
}

void Graph::enum_related_w(vector<int> clique, int delta, ofstream& relative_ofs){
    int* acnt= new int[attr_size];
		if(debug_is_on){
			for(int i=0; i<clique.size(); i++)
				printf("%d ", clique[i]);
			//printf("\n");
		}
	if(clique.size()<threshold * attr_size) return;
	for(int i=0; i<attr_size; i++) acnt[i] = 0;
    for(int i=0; i< clique.size(); i++){
        acnt[attribute[clique[i]]]++;
    }
    int attr_min= acnt[0], attr_max;
    for(int i=1;i <attr_size; i++){
        if(attr_min > acnt[i]) attr_min = acnt[i];
    }
    attr_max = attr_min + delta;
    vector<int> partial_clique;
    //find those less than attr_max;
    vector<int> enum_attrs;
    vector<vector<int>> enum_nodes; 
    for(int i=0; i<attr_size; i++){
        if(acnt[i] <= attr_max){
            for(int j=0; j<clique.size(); j++){
                if(attribute[clique[j]] == i) partial_clique.push_back(clique[j]);
            }
        }else{
            enum_attrs.push_back(i);
            vector<int> nodes;
            for(int j=0; j<clique.size(); j++){
                if(attribute[clique[j]] == i) nodes.push_back(clique[j]);
            }
            enum_nodes.push_back(nodes);
        }
    }
    if(partial_clique.size() == clique.size()){
			relative_num ++;
			relative_ofs<<clique.size()<<" ";
		//	for(int i=0; i<clique.size(); i++)
		//		relative_ofs<<clique[i]<<"("<<attribute[clique[i]]<<")"<<" ";
		//	relative_ofs<<"\n";
			printClique(clique, relative_ofs);
    }else{
			if(debug_is_on){
				printf("attr_min=%d, attr_max=%d\n", attr_min, attr_max);
				for(int ii=0; ii<enum_attrs.size(); ii++)
					printf("%d ",enum_attrs[ii]);
				printf("\n");
				printf("partial clique=");
				for(int ii=0; ii<partial_clique.size(); ii++)
					printf("%d ", partial_clique[ii]);
				printf("\n");
			}

        deep_search(partial_clique, enum_attrs, enum_nodes, 0, 0, attr_max, relative_ofs,0);
    }

}
void Graph::FindRelatedFairClique_W(const char* ord, int delta){


		ofstream relative_ofs("data/relative_cliques_weak.txt");
		printf("start finding relative-weak fair cliques\n");
	  Timer t;
    FindClique(ord);
		printf("size of weak fair cliques is %d\n", ResultSet.size());
		relative_num = 0;
    assert(ResultSet.size()!=0);
    for(int i=0; i<ResultSet.size(); i++){
				sort(ResultSet[i].begin(), ResultSet[i].end());
				vector<int> clique;
				for(vector<int>::iterator it = ResultSet[i].begin(); it!=ResultSet[i].end(); it++)
					clique.push_back(*it);
        enum_related_w(clique, delta, relative_ofs);
    }
		printf("relative weak num =%d\n", relative_num);
		printf("time for relative fair cliques(using weak): %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
		relative_ofs.close();
}

void Graph::RFair_Backtrack(vector<int> R, vector<int>* candidates, int candidate_siz, vector<int> X, vector<int> att_cnt, int attr_max, int delta, int tar_attr, ofstream& relative_ofs){
    vector<int> newR = R;
    change_mem(sizeof(int)*newR.capacity(), 1);
	if(debug_is_on){
		printf("cur R = ");
		for(int i=0; i<R.size(); i++)
			printf("%d ",R[i]);
		printf("\n");
		printf("X =" );
		for(int i=0; i<X.size(); i++) printf("%d ",X[i]);
		printf("\n");
	}
	if(attr_max == -1 && candidates[tar_attr].size() == 0){
		attr_max = delta + att_cnt[tar_attr];
        for(int i=0 ; i<attr_size; i++){
            if(att_cnt[i]>=attr_max){
                candidate_siz -= candidates[i].size();
                candidates[i].clear();
            }
        }
		if(debug_is_on) printf("attr_max become %d candidate size=%d\n", attr_max,candidate_siz);
	}
	if(candidate_siz == 0){
        assert(attr_max > 0);
		int toPrint=true;
		if(X.size() != 0){
			//如果加一个进去，却不会破坏它的公平性，就不是答案
			int attr_min = attr_max - delta;
			for(int i=0; i<X.size(); i++){
				if(att_cnt[attribute[X[i]]]==attr_min ||att_cnt[attribute[X[i]]]+1 <= attr_max){
					toPrint = false;break;
				}
			}
		}
		if(toPrint){
			if(debug_is_on) printf("print\n");
            if(R.size()>=threshold*attr_size){
                relative_num ++;
                relative_ofs<<R.size()<<" ";
			    printClique(R, relative_ofs);
            }
		}
	}
	int* newCnt=new int[attr_size];
    change_mem(sizeof(int)*attr_size,1);
	if(attr_max != -1){
        bool find = false;
		for(int tar = tar_attr; tar<=tar_attr+attr_size; tar++){
			if(candidates[tar%attr_size].size() == 0) continue;
			if(att_cnt[tar%attr_size]<attr_max){
				tar_attr = tar%attr_size; 
                find = true; break;
			}
		}
        if(!find) return ;

	}

    assert(tar_attr>=0 && tar_attr <attr_size);
	if(candidates[tar_attr].size() != 0){
		att_cnt[tar_attr]++;
    }
    if(debug_is_on)
	    printf("tar_attr= %d, att_cnt[%d]=%d\n", tar_attr, tar_attr, att_cnt[tar_attr]);
	for(int i=0; i<candidates[tar_attr].size(); i++){
		int cur=candidates[tar_attr][i];
		newR.push_back(cur);
		for(int j=offset[cur]; j<pend[cur]; j++){
            nvis[edge_list[j]]=1;
        }
		//construct new candidates
		vector<int> newC[attr_size]; 
		for(int j=0; j<attr_size; j++) newCnt[j]=0;
		candidate_siz = 0; 
		int pre_candidate_siz=0;
		int min_candidate_siz = left.size();
		for(int j=0; j<attr_size; j++){
				if(attr_max != -1 && att_cnt[j] >= attr_max) {
					continue;
				}
				if(j==tar_attr){
						for(int k=i+1; k<candidates[j].size(); k++){
								if(nvis[candidates[j][k]]){
										newC[j].push_back(candidates[j][k]);
										newCnt[j]++;
										candidate_siz++;
								}
						}
				}
				else{
						for(int k=0; k<candidates[j].size(); k++){
								if(nvis[candidates[j][k]]){
										newC[j].push_back(candidates[j][k]);
										newCnt[j]++;
										candidate_siz++;
								}
						}
				}
				if(candidate_siz-pre_candidate_siz >0 && candidate_siz-pre_candidate_siz<min_candidate_siz) 
					min_candidate_siz = candidate_siz- pre_candidate_siz;
				pre_candidate_siz = candidate_siz;
		}
        for(int j =0; j<attr_size; j++){
            change_mem(newC[j].capacity()*sizeof(int), 1);
        }
		int cut_flag=0;
		for(int j=0; j<attr_size; j++){
            if(att_cnt[j]+newCnt[j]<threshold){
                    cut_flag=1;break;
            }
		}
		if(cut_flag||candidate_siz+newR.size()<threshold*attr_size){
		//if(cut_flag||(min_candidate_siz+1)*attr_size+newR.size()<threshold*attr_size){
            for(int j=offset[cur]; j<pend[cur]; j++){
                    nvis[edge_list[j]]=0;
            }
            newR.pop_back();
            for(int j =0; j<attr_size; j++){
                change_mem(newC[j].capacity()*sizeof(int), 0);
            }
            continue;
        }
		vector<int> newX;newX.clear();
        change_mem(sizeof(newX)+sizeof(int)*newX.capacity(), 1);
		for(int j=0; j<X.size(); j++){
				if(nvis[X[j]]) newX.push_back(X[j]);
		}
		/**recover**/
		for(int j=offset[cur]; j<pend[cur]; j++){
				nvis[edge_list[j]]=0;
		}
		/****/
		RFair_Backtrack(newR, newC, candidate_siz, newX, att_cnt, attr_max, delta, (tar_attr+1)%attr_size, relative_ofs);
		X.push_back(cur);
		newR.pop_back();
        change_mem(sizeof(newX)+sizeof(int)*newX.capacity(), 0);
        for(int j =0; j<attr_size; j++){
            change_mem(newC[j].capacity()*sizeof(int), 0);
        }
	}
	delete[] newCnt;
    change_mem(sizeof(int)*attr_size,0);
    change_mem(sizeof(int)*newR.capacity(), 0);
}

void Graph::FindRelatedFairClique_S(const char* ord, int delta){
    ofstream relative_ofs("data/relative_cliques_strong.txt");
    Timer t;
    if(attr_size == 2){
        printf("use relative fair reduction\n");
        Relative_fair_reduction();
    }
    else colorfulcore();
    clear_mem();
    int* vis = new int[n];
    change_mem(sizeof(int)*n, 1);
    nvis=new int[n];
    change_mem(sizeof(int)*n, 1);
    for(int i=0; i<n; i++) nvis[i]=0;
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<left.size(); i++)
        vis[left[i]]=1;
    relative_num = 0;
    vector<int> R, X;
    vector<int> candidates[attr_size];
    vector<int> att_cnt;
    for(int i=0; i<attr_size; i++)
        att_cnt.push_back(0);
    int candidate_siz=0;
    relative_num = 0;
    for(int i=0; i<n; i++){
    if(vis[i]){
        component.clear();
        get_connected_component(i, vis);
        change_mem(sizeof(int)*component.size(), 1);
        R.clear(); X.clear();
        for(int j=0; j<attr_size; j++)
            candidates[j].clear();
        for(int j=0; j<attr_size; j++)
            att_cnt[j] = 0;
        printf("siz=%d\n", component.size());
        if(strcmp(ord,"sort")==0)
                sort(component.begin(), component.end());
        else if(strcmp(ord, "colorful")==0)
                component=GetColorfulOrdering();
        else if(strcmp(ord, "degeneracy") ==0)
                component = degeneracy_ordering(component);
        else if(strcmp(ord, "relative") == 0)
                component = GetRelativeOrdering(component);
        else if(strcmp(ord, "anonymous")==0)
                component=GetFairnessOrdering(attr_size);
        for(int j=0; j<component.size(); j++){
            candidates[attribute[component[j]]].push_back(component[j]);
            candidate_siz++;
        }
        printf("start search!\n");
        RFair_Backtrack(R, candidates, candidate_siz, X, att_cnt, -1, delta, 0, relative_ofs);
        change_mem(sizeof(int)*component.size(), 0);
        }
    }
    printf("max mem =%d\n", max_mem);
    printf("relative weak num =%d\n", relative_num);
    printf("time for relative fair cliques(using strong): %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::Relative_fair_reduction(){
    colorfulcore();  //k->k/
    clear_mem();
    printf("left size1=%d\n",left.size());
    Timer tt;
    int* vis = new int[n];
    change_mem(sizeof(int)*n, 1);
    for(int i=0; i<n; i++) vis[i]=1;
    for(int i=0; i<left.size(); i++) vis[left[i]]=0;

    fairness_d=new int[n];
    change_mem(sizeof(int)*n, 1);
    for(int i=0; i<n; i++) fairness_d[i]=0;
    //calculate groups
    int** cntGroup=new int*[n];
    change_mem(sizeof(int)*n*3, 1);
    for(int i=0; i<n; i++){
        cntGroup[i]= new int[3];
        for(int j=0; j<3; j++){
            cntGroup[i][j]=0;
        }
    }

    //c1-> number of groups of color 0,
    //c2-> number of groups of color 1;
    int c1, c2, cboth, ccnt, typ;

    for(int i=0; i<n; i++){
        if(vis[i]) continue;
        c1=c2=cboth=0;
        int num=0;
        for(int j=0; j<=max_color; j++){
            for(int k=0; k<attr_size; k++)
                cntGroup[j][k]=0;
        }
        for(int j=offset[i]; j<pend[i]; j++){
            if(!vis[edge_list[j]])
                ++cntGroup[color[edge_list[j]]][attribute[edge_list[j]]];
        }
        for(int j=0; j<max_color; j++){
            ccnt=0;typ=0;
            for(int k=0; k<attr_size; k++){
                if(cntGroup[j][k]){
                    ccnt++;
                    typ=k;
                    num=cntGroup[j][k];
                }
            }
            if(ccnt==1){
                if(typ==0){ 
                    c1++;
                }
                else{
                    c2++;
                }
            }
            else if(ccnt>1){
                cboth++;
            }
        }
        if(attribute[i]==0){
            c1++;
        }
        else c2++;

        //calculate fairness degree
        if(c1 < threshold) cboth -= threshold - c1;
        if(c2 < threshold) cboth -= threshold -c2;
        fairness_d[i] = cboth;
    }

    //queue<int> Q;
    printf("threshold=%d\n", threshold);
    int deletenode=0;
    int updatenode=0;
    for(int i=0; i<n; i++){
        if(!vis[i]&&fairness_d[i]<0){
            updatenode++;
            vis[i]=1;
        }
    }
    while(updatenode){
        updatenode=0;
        for(int i=0; i<n; i++){
            if(vis[i]) continue;
            c1=c2=cboth=0;
            int num=0;
            for(int j=0; j<=max_color; j++){
                for(int k=0; k<attr_size; k++)
                    cntGroup[j][k]=0;
            }
            for(int j=offset[i]; j<pend[i]; j++){
                if(!vis[edge_list[j]])
                    ++cntGroup[color[edge_list[j]]][attribute[edge_list[j]]];
            }
            for(int j=0; j<max_color; j++){
                ccnt=0;typ=0;
                for(int k=0; k<attr_size; k++){
                    if(cntGroup[j][k]){
                        ccnt++;
                        typ=k;
                        num=cntGroup[j][k];
                    }
                }
                if(ccnt==1){
                    if(typ==0){ 
                        c1++;
                    }
                    else{
                        c2++;
                    }
                }
                else if(ccnt>1){
                    cboth++;
                }
            }
            if(attribute[i]==0){
                c1++;
            }
            else c2++;

            //calculate fairness degree
            if(c1 < threshold) cboth -= threshold - c1;
            if(c2 < threshold) cboth -= threshold -c2;
            fairness_d[i] = cboth;
            //printf("fairness degree[%d]=%d\n", i, fairness_d[i]);
            if(fairness_d[i]<0){
                updatenode++;
                vis[i]=1;
            }
        }
    }
    //construct subgraph
    construct_sub_graph(vis);
    delete[] vis;
    change_mem(sizeof(int), 0);
    delete[] cntGroup;
    change_mem(sizeof(int)*n*3, 0);
    printf("max mem for strong reduction=%lld, left node =%d\n", max_mem, left.size());
    printf("time for second reduction: %s (microseconds)\n", Utility::integer_to_string(tt.elapsed()).c_str());
}

vector<int> Graph::GetRelativeOrdering(vector<int> array){
	printf("get relative order\n");
    int* vis = new int[n];
    for(int i=0; i<n; i++) vis[i]=1;
    for(int i=0; i<component.size(); i++) vis[component[i]]=0;

    int** cntGroup=new int*[n];
    for(int i=0; i<n; i++){
        cntGroup[i]= new int[3];
        for(int j=0; j<3; j++){
            cntGroup[i][j]=0;
        }
    }
    int* queue_n= new int[component.size()];
    for(int i=0; i<component.size(); i++) queue_n[i]=i;
    int* queue_fair= new int[component.size()];
    for(int i=0; i<component.size(); i++) queue_fair[i]=fairness_d[component[i]];
    int* idx= new int[n];
    for(int i=0; i<component.size(); i++)
        idx[component[i]]=i;

    vector<int> peeling_order; peeling_order.clear();
    ListLinearHeap *heap = new ListLinearHeap(component.size(), component.size());//nodes, max degree
    heap->init(component.size(), component.size(), queue_n, queue_fair);//number of nodes, max degree, start position of arr, start pos of degreesa
    change_local_mem(heap->get_max_occupy(), 1);
    for(ui i = 0;i < component.size();i ++) {
        ui u, key;
        heap->pop_min(u, key);
        int node=component[u];
        peeling_order.push_back(node);
        vis[node] = 1;
        for(ui nn = offset[node];nn < pend[node];nn ++) if(vis[edge_list[nn]] == 0) {
            int nei=edge_list[nn];
        //c1-> number of groups of color 0,
        //c2-> number of groups of color 1;
            int old_fairness=fairness_d[nei];
            int c1, c2, cboth, ccnt, typ;
            c1=c2=cboth=0;
            int num=0;
            for(int j=0; j<=max_color; j++){
                for(int k=0; k<attr_size; k++)
                    cntGroup[j][k]=0;
            }
            for(int j=offset[nei]; j<pend[nei]; j++){
                if(!vis[edge_list[j]])
                    ++cntGroup[color[edge_list[j]]][attribute[edge_list[j]]];
            }
            for(int j=0; j<max_color; j++){
                ccnt=0;typ=0;
                for(int k=0; k<attr_size; k++){
                    if(cntGroup[j][k]){
                        ccnt++;
                        typ=k;
                        num=cntGroup[j][k];
                    }
                }
                if(ccnt==1){
                    if(typ==0) c1++;
                    else c2++;
                }
                else if(ccnt>1){
                    cboth++;
                }
            }
            if(attribute[nei]==0) c1++;
            else c2++;
            //calculate fairness degree

            if(c1 < threshold) cboth -= threshold - c1;
            if(c2 < threshold) cboth -= threshold -c2;
            fairness_d[nei] = cboth;
            if(fairness_d[nei] >= 0)
                heap->decrement(idx[nei], old_fairness-fairness_d[nei]);//decrease the degree by 1
        }
    }
    delete heap;
    delete[] cntGroup;
    delete[] vis;
    delete[] queue_n;
    delete[] queue_fair;
    delete[] idx;
    return peeling_order;
}
