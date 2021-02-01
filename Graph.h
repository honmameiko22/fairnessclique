#ifndef GRAPH_H_
#define GRAPH_H_
#include<bits/stdc++.h>

using namespace std;

class Graph{
    private:
    int n, m;//num of nodes and edges
    int attr_size;//number of attributes
    int threshold;//given
    int max_color;
    string dir;
    int* offset;//index of nodes in edge_list
    int* edge_list;//store all neighbors
    int* attribute;//the attribute of node
    int* pend;    
    int* nei_vis;
    int* fairness_d;
    long long now_mem=0;
    long long max_mem=0;
    long long max_local_mem=0;
    long long now_local_mem=0;
    vector<int> left;
    vector<int> component;
    vector<int> union_X;
    vector< vector<int> > ResultSet;
    
    int* peeling_idx;
    int* color;

    int* X_vis;
    int* nvis;
    int* rvis;


    public:
    Graph(const char* _dir);//init some values
    ~Graph();
    void ReadGraph(const char* inputname, const char* attr, int, const char* mode);
    void colorfulcore();
    void FindClique(const char*);
    void FindStrongClique(const char*);
    void Strong_reduction_basic();
    void Strong_reduction();
    void Strong_reduction_Heuristic();
    void Strong_reduction_Three();
    void VerifyWClique(const char* inputname);
    void VerifySClique(const char* inputname);
    vector<int> core_decomposition(int);
    vector<int> core_decomposition_cut();
    vector<int> GetFairnessOrdering(int);
    vector<int> GetDegreeOrdering();
    vector<int> FairnessDegeranacyOrdering();
    vector<int> GetColorfulFairnessOrdering();
    vector<int> GetColorfulFairnessOrdering_Heuristic();
    vector<int> GetColorfulFairnessOrdering_Heuristic_2();
    void Backtrack_Strong(vector<int>&, vector<int>* , vector<int>&, int);
    void BBE(vector<int>&, vector<int>&);
    void Backtrack_Strong_modify(vector<int>&, vector<int>* , vector<int>&, int);
    void Backtrack_Weak(vector<int>&, vector<int>&, vector<int>&, int*, int);
    void Backtrack_Weak_improve(vector<int>&, vector<int>&, vector<int>&, int*, int);
    long long get_max();

    private:
    void change_mem(long long tmp, int flag){
        if(flag){
            now_mem+=tmp;
            max_mem = now_mem>max_mem? now_mem: max_mem;
        }
        else now_mem-=tmp;
    }
    void clear_mem(){
        now_mem=0;
        max_mem=0;
    }
    void change_local_mem(long long tmp, int flag){
        if(flag){
            now_local_mem+=tmp;
            max_local_mem=now_local_mem>max_local_mem?now_local_mem:max_local_mem;
        }
        else now_local_mem-=tmp;
    }
    void clear_local_mem(){
        now_local_mem=0;
    }
    int Verify_largest(vector<int>& X){
        vector<int> set_[attr_size];
        for(int i=0; i<attr_size; i++) set_[i].clear();

        for(int j=0; j<X.size(); j++){
            set_[attribute[X[j]]].push_back(X[j]);
        }
        long long tmp_sum=0;
        for(int j=0; j<attr_size; j++) 
            tmp_sum+=sizeof(int)*(set_[j].capacity());
        change_mem(sizeof(int)*tmp_sum, 1);
        int empty=0;
        for(int i=0; i<attr_size; i++){
            if(set_[i].size()==0) empty++;
        }
        int have=1;//i->can be added, 0->cannot be added
        if(empty==0){
            vector<vector<int> >record, swap_record;
            record.clear();
            for(int i=0; i<set_[0].size(); i++){
                vector<int> vec{set_[0][i]};
                record.push_back(vec);
            }
            for(int level=1; level<attr_size; level++){
                swap_record.clear();
                for(int i=0; i<set_[level].size(); i++){
                    int nod=set_[level][i];
                    for(int j=offset[nod]; j<pend[nod]; j++)
                        X_vis[edge_list[j]]=1;
                    for(int j=0; j<record.size(); j++){
                        int allin=1;
                        for(int k=0; k<record[j].size(); k++){
                            if(X_vis[record[j][k]]==0){
                                allin=0;break;
                            }
                        }
                        if(allin){//add to clique (record)
                            vector<int> temp=record[j];
                            temp.push_back(nod);
                            swap_record.push_back(temp);
                        }
                    }
                    for(int j=offset[nod]; j<pend[nod]; j++)
                        X_vis[edge_list[j]]=0;
                }
                record=swap_record;
                if(record.empty()) break;
            }
            if(!record.empty()) have=0;
            change_mem(sizeof(int)*record.capacity()+sizeof(record),1);
            change_mem(sizeof(int)*swap_record.capacity()+sizeof(swap_record),1);
            
            change_mem(sizeof(int)*record.capacity()+sizeof(record),0);
            change_mem(sizeof(int)*swap_record.capacity()+sizeof(swap_record),0);
        }

        change_mem(sizeof(int)*tmp_sum, 0);
        return have;
    }

    void construct_sub_graph(int* vis){
        int startpos=0;
        left.clear();
        for(int i=0; i<n; i++){
            if(!vis[i]){
                left.push_back(i);
                int offset_start=startpos;
                for(int j=offset[i]; j<pend[i]; j++){
                    if(!vis[edge_list[j]]){
                        edge_list[startpos++]=edge_list[j];
                    }
                }
                offset[i]=offset_start;
                pend[i]=startpos;
            }
            else{
                offset[i]=pend[i]=0;
            }
        }
        offset[n]=pend[n]=0;
        //printf("2--< the left nodes of the graph is %d > \n", left.size());
        //printf("2--< the left edges of the graph is %d >\n", startpos);
    }


    void get_connected_component(int root, int* vis){
        stack<int> Q;
        while(!Q.empty()) Q.pop();
        Q.push(root);
        vis[root]=0;
        int maxQ=0;
        while(!Q.empty()){
            if(Q.size()>maxQ) maxQ=Q.size();
            int cur =Q.top(); Q.pop();
            component.push_back(cur);
            for(int i=offset[cur]; i<pend[cur]; i++){
                if(vis[edge_list[i]]){
                    Q.push(edge_list[i]);
                    vis[edge_list[i]]=0;
                }
            }
        }
        change_mem(sizeof(component)+ sizeof(int)*component.capacity(), 1);
        change_mem(sizeof(Q)+sizeof(int)*maxQ, 1);
        change_mem(sizeof(Q)+sizeof(int)*maxQ, 0);
    }
};


#endif