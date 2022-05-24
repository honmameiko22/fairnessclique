#include "Utility.h"
#include "Graph.h"
#include "LinearHeap.h"
#include "Timer.h"
using namespace std;

Graph::Graph(const char* _dir){//init some values
    dir = string(_dir);
    n = 0;
    m = 0;
    threshold = 0;
    attr_size=0;
    max_color=0;
    offset = nullptr;
    edge_list = nullptr;
    attribute = nullptr;
    peeling_idx = nullptr;
    color = nullptr;
    pend = nullptr;
    nei_vis= nullptr;
    nvis = nullptr;
    rvis = nullptr;
    X_vis = nullptr;
    peeling_idx= nullptr;
    fairness_d= nullptr;
    left.clear();
    component.clear();
    union_X.clear();
    ResultSet.clear();
}
Graph::~Graph(){
    if(offset != nullptr){
        delete[] offset;
        offset = nullptr;
    }
    if(edge_list != nullptr){
        delete[] edge_list;
        edge_list = nullptr; 
    }
    if(attribute != nullptr){
        delete[] attribute; 
        attribute = nullptr;
    }
    if(pend != nullptr){
        delete[] pend;
        pend = nullptr;
    }
    if(peeling_idx != nullptr){
        delete[] peeling_idx;
        peeling_idx=nullptr;
    }
    if(color != nullptr){
        delete[] color;
        color = nullptr;
    }
    if(X_vis!=nullptr){
        delete[] X_vis;
        X_vis=nullptr;
    }
    if(nvis != nullptr){
        delete[] nvis;
        nvis = nullptr;
    }
    if(rvis != nullptr){
        delete[] rvis;
        rvis = nullptr;
    }
}
void Graph::ReadGraph(const char* inputname, const char* attr, int th, const char* mode){
    threshold=th;
    printf("Start reading graph from %s\n", inputname);
    FILE *f = fopen(inputname, mode);
    if(f == nullptr){
        printf("Can not open the graph file!\n");
        exit(1);
    }
    vector< pair<int, int> > pairs;
    int from, to;
    while(fscanf(f, "%d%d", &from, &to) == 2){//read undirected edges, ordered edges
        if(from == to) continue;
        pairs.push_back(make_pair(from, to));
    }
    for(int i = 0; i < pairs.size(); i++){
        if(pairs[i].first > n) n = pairs[i].first;
        if(pairs[i].second > n) n = pairs[i].second;
        if(i > 0 && pairs[i].first == pairs[i-1].first && pairs[i].second == pairs[i-1].second) printf("Exist self loop!\n");  
    }
    n++;//range from [0, n-1];
    m=pairs.size()+1;
    change_mem(sizeof(pairs)+2*sizeof(int)*pairs.size(), 1);
    printf("< Number of nodes= %d, number of edges= %d >\n", n, m);

    //store edges
    edge_list = new int[m+1];
    change_mem((m+1)*sizeof(int), 1);
    offset = new int[n+1];
    change_mem((n+1)*sizeof(int), 1);
    offset[0] = 0;
    int indx = 0;
    for(int i = 0; i < n ; i++){
        offset[i+1] = offset[i];
        while(indx < pairs.size() && pairs[indx].first == i) edge_list[offset[i+1]++] = pairs[indx++].second;
    }
    fclose(f);
    //read attribute
    //printf("Start reading attributes from %s\n", attr);
    attribute=new int[n];
    change_mem(n*sizeof(int),1);
    f = fopen(attr, mode);
    if(f == nullptr){
        printf("Can not open the attribute file\n");
        exit(1);
    }
    int node, a;
    attr_size=0;
    while(fscanf(f, "%d%d", &node, &a)==2){
        if(node>=n){
           printf("attribute error! node=%d\n",node);
           exit(1);
        }
        attribute[node]=a;
        if(a>attr_size) attr_size=a;
    }
    attr_size++;
    fclose(f);
    change_mem(sizeof(pairs)+2*sizeof(int)*pairs.size(), 0);
    printf("mem of graph = %lld\n", (long long)((m+1)*sizeof(int)+(n+1)*sizeof(int)));
    printf("mem of attribute=%ld\n", (long long)n*sizeof(int));
    printf("read done\n");
}

void Graph::colorfulcore(){
    //compute colorful degree
    Timer t;
    printf("start colorfulcore...\n");
    int* cvis=new int[n];
    change_mem(sizeof(int) * n, 1);
    color = new int[n];
    change_mem(sizeof(int) * n, 1);
    int* head = new int[n];
    change_mem(sizeof(int) * n, 1);
    int* nxt = new int[n];
    change_mem(sizeof(int) * n, 1);
    int max_degree=0;
    int* degree = new int[n];
    change_mem(sizeof(int) * n, 1);
    pend = new int[n+1];
    change_mem(sizeof(int) * (n+1), 1);
    int* vis = new int [n];//vis==1->deleted
    change_mem(sizeof(int) * n, 1);
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<n+1; i++) pend[i]=0;


    for(int i=0; i<n; i++){
        degree[i] = offset[i+1]-offset[i];
        pend[i]=offset[i+1];
        head[i] = n;
    }
    int cut_thre=threshold*attr_size-1;
    int first_reduction=0;
    for(int i=0; i<n; i++){
        if(degree[i]<cut_thre){
            vis[i]=1;first_reduction++;
            continue;
        }
        nxt[i]=head[degree[i]];
		head[degree[i]]=i;
        if(degree[i] > max_degree) max_degree=degree[i];
    }
    delete[] degree;
    change_mem(sizeof(int) * n, 0);
    //printf("first reduction=%d\n", first_reduction);
    
    //printf("Coloring...\n");
    for(int i=0; i<n; i++) cvis[i]=0;
    for(int i = 0;i < n;i ++) color[i] = n;//n means illegal
	//decide the color of vertices
	max_color = 0;
	for(int ii=max_degree; ii>=1; ii--){
		for(int jj=head[ii]; jj!=n; jj=nxt[jj]){
			int u = jj;
			for(int j = offset[u];j < offset[u+1];j ++) {
                int c = color[edge_list[j]];
                if(c != n) {
                    cvis[c] = 1;
                }
			}

			for(int j = 0;;j ++){
                if(!cvis[j]) {
                    color[u] = j;
                    if(j > max_color) max_color = j;
                    break;
			    }
            }
			for(int j = offset[u];j < offset[u+1];j ++) {
                int c = color[edge_list[j]];
                if(c != n) cvis[c] = 0;
			}
		}
	}
    max_color++;
    delete[] cvis;
    change_mem(sizeof(int) * n, 0);
    delete[] head;
    change_mem(sizeof(int) * n, 0);
    delete[] nxt;
    change_mem(sizeof(int) * n, 0);
    #ifdef DEBUG
    for(int i=0; i<n; i++){
        printf("color[%d]=%d\n", i, color[i]);
    }
    #endif
    //compute d_attribute[node]

    //printf("max color=%d\n", max_color);
    int** colorful_d = new int*[n];
    int*** colorful_r =  new int **[n];
    change_mem(sizeof(int) * n*attr_size, 1);
    change_mem(sizeof(int) * n*attr_size*max_color, 1);
    for(int i=0; i<n; i++){
        colorful_d[i] = new int[attr_size];
        for(int j=0; j<attr_size; j++)
            colorful_d[i][j]=0;
    }
    for(int i=0; i<n; i++){
        colorful_r[i] = new int*[attr_size];
        for(int j=0; j<attr_size; j++){
            colorful_r[i][j] = new int[max_color];
            for(int k=0; k<max_color; k++)
                colorful_r[i][j][k]=0;
        }
    }
    for(int i=0; i<n; i++){
        for(int j=offset[i]; j<offset[i+1]; j++){
            int neighbor=edge_list[j];
            if(vis[neighbor]==1) continue;
            if((colorful_r[i][attribute[neighbor]][color[neighbor]]++)==0){
                colorful_d[i][attribute[neighbor]]++;
            }
        }
    }
    printf("compute colorful degree done!\n");

    //compute colorful core
    queue<int> Q;//Q is the set to delete
    int deletenode=0;
    for(int i=0; i<n; i++){
        deletenode=0;
        if(vis[i]==1) continue;
        if(colorful_d[i][attribute[i]]<threshold-1)
            deletenode=1;
        if(deletenode==0){
            for(int j=0; j<attribute[i]; j++){
                if(colorful_d[i][j]<threshold){
                    deletenode=1; break;
                }
            }
            if(deletenode==0){
                for(int j=attribute[i]+1; j<attr_size; j++){
                    if(colorful_d[i][j]<threshold){
                        deletenode=1; break;
                    }
                }
            }
        }
        if(deletenode == 1){
            Q.push(i);
            vis[i]=1;
        }
    }
    int maxQ=Q.size();
    while(!Q.empty()){
        if(Q.size()>maxQ) maxQ=Q.size();
        int cur=Q.front();
        int mycolor=color[cur];
        Q.pop(); 
        for(int i=offset[cur]; i<offset[cur+1]; i++){
            int neighbor=edge_list[i];
            if(!vis[neighbor]){
                if(--colorful_r[neighbor][attribute[cur]][color[cur]] <= 0){
                    colorful_d[neighbor][attribute[cur]]--;
                }
                deletenode=0;

                if(colorful_d[neighbor][attribute[neighbor]]<threshold-1)
                    deletenode=1;
                if(deletenode==0){
                    for(int j=0; j<attribute[neighbor]; j++){
                        if(colorful_d[neighbor][j]<threshold){
                            deletenode=1; break;
                        }
                    }
                    if(deletenode==0){
                        for(int j=attribute[neighbor]+1; j<attr_size; j++){
                            if(colorful_d[neighbor][j]<threshold){
                                deletenode=1; break;
                            }
                        }
                    }
                }

                if(deletenode==1){
                    Q.push(neighbor);
                    vis[neighbor]=1;//other wise a node will be push multiple times
                }
            }
        }
    }
    change_mem(sizeof(Q)+sizeof(int) * maxQ, 1);
    printf("max mem for colorful core=%lld\n",
    ((long long)sizeof(int))*(3*n+1+n*attr_size+n*attr_size*max_color+maxQ)+sizeof(Q)
    );
    left.clear();
    delete[] colorful_r;
    delete[] colorful_d;
    change_mem(sizeof(int) * n*attr_size, 0);
    change_mem(sizeof(int) * n*attr_size*max_color, 0);
    int lastp;
    
    for(int i=0; i<n; i++){
        if(!vis[i]){
            left.push_back(i);
        }
    }
    change_mem(sizeof(left)+ sizeof(int)*left.capacity(), 1);
    printf("< the left nodes of the graph is %d > \n", left.size());

    
    //reconstruct graph
    int startpos=0;
    for(int i=0; i<n; i++){
        if(!vis[i]){
            int offset_start=startpos;
            for(int j=offset[i]; j<offset[i+1]; j++){
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
    printf("< the left edges of the graph is %d >\n", startpos);
    printf("second reduction=%d\n", n-first_reduction-left.size());
    #ifdef DEBUG
    for(int i=0; i<n; i++){
        printf("i=%d: ",i);
        for(int j=offset[i]; j<pend[i]; j++){
            printf(" %d", edge_list[j]);
        }
        printf("\n");
    }
  #endif
    //use pend instead of offset[i+1]
    delete[] vis;
    change_mem(sizeof(int)*n, 0);
    change_mem(sizeof(Q)+ sizeof(int)*maxQ, 0);
    printf("time for colorful-core: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
}
int acnt=0;

void Graph::FindClique(const char* ord){
    Timer t;
    //if(attr_size == 2){
    //    printf("use relative fair reduction\n");
    //    Relative_fair_reduction();
    //}
    //else colorfulcore();
    colorfulcore();
    clear_mem();
    //printf("colorfulcore done!\n");
    printf("time for colorfulcore: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
    vector<int> R, X;
    change_mem(sizeof(R), 1);
    change_mem(sizeof(X), 1);
    nvis=new int[n];
    change_mem(sizeof(int)*n ,1);
    rvis=new int[n];
    change_mem(sizeof(int)*n ,1);
    int* cntC=new int[attr_size];//cntC is the count of node of a centain attribute in candidate set
    change_mem(sizeof(int)*attr_size ,1);
    int* vis = new int[n];//vis denotes the set of nodes unsearched
    change_mem(sizeof(int)*n ,1);
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<n; i++) nvis[i]=0;
    for(int i=0; i<n; i++) rvis[i]=0;
    //change into smaller parts
    for(int i=0; i<left.size(); i++) 
        vis[left[i]]=1;
    ResultSet.clear();
    long long max_mem_for_ordering=0;
    
    for(int i=0; i<n; i++){
        if(vis[i]){
            component.clear();
            get_connected_component(i, vis);
            printf("component size=%d\n", component.size());
            #ifdef DEBUG
            for(int i=0; i<component.size(); i++){
                printf("%d ", component[i]);
            }
            printf("\n");
            getchar();
            #endif
            int csize=component.size();
            for(int j=0; j<attr_size; j++) cntC[j]=0;
            for(int j=0; j<csize; j++){//in the beginning, all nodes belong to the candidate set
                cntC[attribute[component[j]]]++;
            }
            R.clear();X.clear();
            #ifdef DEBUG
            for(int j=0; j<attr_size; j++){
                printf("cnt[%d]=%d\n", j, cntC[j]);
            }
            #endif
            if(strcmp(ord,"sorted")==0) 
                sort(component.begin(), component.end());
            else if(strcmp(ord,"colorful")==0)   
                component=GetColorfulOrdering();
            else if(strcmp(ord,"degree")==0)
                component=GetDegreeOrdering();
            else if(strcmp(ord, "degeneracy") == 0)
                component = degeneracy_ordering(component);
            else if(strcmp(ord, "anonymous") == 0)
                component=GetFairnessOrdering(attr_size);
            else if(strcmp(ord, "relative") == 0)
                component = GetRelativeOrdering(component);
            long long tmp=(4*n+component.size()*attr_size+component.size()*attr_size*max_color+component.size()*3)*(long long)sizeof(int);
            if(tmp>max_mem_for_ordering) max_mem_for_ordering=tmp;
            Backtrack_Weak(R, component, X, cntC, 1);
        }
    }
    printf("acnt=%d\n", acnt);
    printf("max_mem_for_ordering=%lld\n", max_mem_for_ordering);
    printf("search done!\n");
    printf("tot time: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
    printf("size of result = %d\n", ResultSet.size());
    /*
    ofstream ofs("fair_clique.txt");
    printf("size of result = %d\n", ResultSet.size());
    for(int i=0; i<ResultSet.size(); i++){
        ofs<<ResultSet[i].size()<<" ";
        for(int j=0; j<ResultSet[i].size(); j++){
            ofs<<ResultSet[i][j]<<" ";
        }
        ofs<<"\n";
    }
    ofs.close();
    */
    delete[] cntC;
    delete[] nvis;
    delete[] rvis;
    delete[] vis;
    change_mem(sizeof(int)*n*4 ,0);
}

void Graph::Backtrack_Weak(vector<int>& R, vector<int>& C, vector<int>& X, int* cntC, int depth){//R->current clique, C->candidates ordered, X->enumerated
    
    if(C.size()==0&&X.size()==0&&R.size()>=threshold*attr_size){
        ResultSet.push_back(R);
        return;
    }
    else if(C.size()==0){
        return;
    }
    vector<int> SwapC; SwapC.clear();
    for(int i=0; i<R.size(); i++) rvis[R[i]]=1;
    for(int i=0; i<C.size(); i++){
        if(!rvis[C[i]]) SwapC.push_back(C[i]);
    }
    C=SwapC;
    change_mem(sizeof(SwapC)+ sizeof(int)*SwapC.capacity(), 1);
    for(int i=0; i<R.size(); i++) rvis[R[i]]=0;
    /**decide the order of C**/
    vector<int> newR; newR.clear();
    for(int j=0; j<R.size(); j++)
        newR.push_back(R[j]);
    change_mem(sizeof(newR)+ sizeof(int)*newR.capacity(), 1);
    int* newCnt=new int[attr_size];
    change_mem(sizeof(int)*attr_size, 1);
    for(int i=0; i<C.size(); i++){
        int cur=C[i];//add node"cur" to the current clique
        newR.push_back(cur);
        //construct new candidate set
        vector<int> newC; newC.clear();
        //present the nodes that has been deleted
        for(int j=offset[cur]; j<pend[cur]; j++){
            nvis[edge_list[j]]=1;
        }
        for(int j=0; j<attr_size; j++) newCnt[j]=0;
        for(int j=i+1; j<C.size(); j++){//avoid finding a clique multiple times
            if(nvis[C[j]]){
                newC.push_back(C[j]);
                newCnt[attribute[C[j]]]++;
            }
        }
        change_mem(sizeof(newC)+sizeof(int)*newC.capacity(), 1);
        int cut_flag=0;
        int* cnt_current=new int[attr_size];
        change_mem(sizeof(int)* attr_size, 1);
        for(int j=0; j<attr_size; j++)
            cnt_current[j]=0;
        for(int j=0; j<newR.size(); j++)
            cnt_current[attribute[newR[j]]]++;

        for(int j=0; j<attr_size; j++){
            if(cnt_current[j]+newCnt[j]<threshold){
                cut_flag=1;break;
            }
        }

        change_mem(sizeof(int)* attr_size, 0);
        delete[] cnt_current;

        if(cut_flag||newC.size()+newR.size()<threshold*attr_size){
            for(int j=offset[cur]; j<pend[cur]; j++){
                nvis[edge_list[j]]=0;
            }
            newR.pop_back();

            change_mem(sizeof(newC)+sizeof(int)*newC.capacity(), 0);
            continue;
        }

        vector<int> newX;newX.clear();
        for(int j=0; j<X.size(); j++){
            if(nvis[X[j]]) newX.push_back(X[j]);
        }
        change_mem(sizeof(newX)+sizeof(int)*newX.capacity(), 1);
        /**recover**/
        for(int j=offset[cur]; j<pend[cur]; j++){
            nvis[edge_list[j]]=0;
        }
        /****/
        Backtrack_Weak(newR, newC, newX, newCnt, depth+1);
        X.push_back(cur);
        newR.pop_back();

        change_mem(sizeof(newC)+sizeof(int)*newC.capacity(), 0);
        change_mem(sizeof(newX)+sizeof(int)*newX.capacity(), 0);
    }
    delete[] newCnt;
    change_mem(sizeof(int)*attr_size, 0);
    change_mem(sizeof(SwapC)+sizeof(int)*SwapC.capacity(), 0);
    change_mem(sizeof(newR)+sizeof(int)*newR.capacity(), 0);
}


void Graph::FindStrongClique(const char* ord){
    Timer t;
    printf("order=%s\n", ord);
    //colorfulcore();
    //Strong_reduction_Three();
    Strong_reduction_basic();
    //Strong_reduction_Heuristic();
    clear_mem();
    nvis=new int[n];change_mem(sizeof(int)*n, 1);
    rvis=new int[n];change_mem(sizeof(int)*n, 1);
    X_vis=new int[n];change_mem(sizeof(int)*n, 1);
    int* vis = new int[n];change_mem(sizeof(int)*n, 1);
    nei_vis=new int[n];change_mem(sizeof(int)*n, 1);
    vector<int> R, X;
    vector<int> candidates[attr_size];
    for(int i=0; i<n; i++) X_vis[i]=0;
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<n; i++) nvis[i]=0;
    for(int i=0; i<n; i++) rvis[i]=0;
    for(int i=0; i<n; i++) nei_vis[i]=0;
    for(int i=0; i<left.size(); i++)
        vis[left[i]]=1;
    ResultSet.clear();
    printf("left size of graph is %d\n", left.size());
    int* head_pos[attr_size];
    change_mem(sizeof(int)*attr_size, 1);
    for(int i=0; i<n; i++){
        if(vis[i]){
            component.clear();
            get_connected_component(i, vis);
            R.clear(); X.clear();
            for(int j=0; j<attr_size; j++)
                candidates[j].clear();
            
            if(component.size()<threshold*attr_size) continue;
            printf("siz=%d\n", component.size());
            if(strcmp(ord,"sort")==0)
                sort(component.begin(), component.end());
            else if(strcmp(ord, "heuristic")==0)
                component=GetColorfulFairnessOrdering_Heuristic_2();
            else if(strcmp(ord, "fairness")==0)
                component=GetColorfulFairnessOrdering();
            else if(strcmp(ord, "anonymous")==0)
                component=GetFairnessOrdering(attr_size);

            //component=GetColorfulFairnessOrdering_Heuristic();
            printf("done\n");
            if(strcmp(ord, "core")!=0)
                for(int j=0; j<component.size(); j++){
                    candidates[attribute[component[j]]].push_back(component[j]);
                }
            
            Backtrack_Strong(R, candidates, X, 0);
        }
    }
    printf("max mem for ordering=%lld\n", max_local_mem);
    printf("size of result = %d\n", ResultSet.size());
   /*
    ofstream ofs("strong_fair_clique.txt");
    printf("size of result = %d\n", ResultSet.size());
    for(int i=0; i<ResultSet.size(); i++){
        ofs<<ResultSet[i].size()<<" ";
        for(int j=0; j<ResultSet[i].size(); j++){
            ofs<<ResultSet[i][j]<<",";
        }
        ofs<<"\n";
    }
    ofs.close();
    
    ofs.open("strong_fair_clique_ordered.txt", std::ofstream::out);
    for(int i=0; i<ResultSet.size(); i++){
        ofs<<ResultSet[i].size()<<" ";
        sort(ResultSet[i].begin(), ResultSet[i].end());
        for(int j=0; j<ResultSet[i].size(); j++){
            ofs<<ResultSet[i][j]<<",";
        }
        ofs<<"\n";
    }
    ofs.close();
    */
    printf("tot time: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
    delete[] rvis;
    delete[] nvis;
    delete[] vis;
    delete[] X_vis;
}
void Graph::Backtrack_Strong(vector<int>& R, vector<int> *candidates, vector<int>& X, int att_node){
    int zero_cnt=0; 
    if(R.size()%attr_size==0&&R.size()>=threshold*attr_size){//check
      //check if 
        change_mem(sizeof(int)*union_X.capacity()+sizeof(union_X), 0);
        union_X.clear();
        for(int j=offset[R[0]]; j<pend[R[0]]; j++){
            nei_vis[edge_list[j]]=1;
        }
        for(int i=1; i<R.size(); i++){
            for(int j=offset[R[i]]; j<pend[R[i]]; j++){
                if(nei_vis[edge_list[j]])
                    nei_vis[edge_list[j]]++;
            }
        }
        for(int j=offset[R[0]]; j<pend[R[0]]; j++){
            if(nei_vis[edge_list[j]]==R.size()) union_X.push_back(edge_list[j]);
            nei_vis[edge_list[j]]=0;
        }
        change_mem(sizeof(int)*union_X.capacity()+sizeof(union_X), 1);
        if(union_X.size()<attr_size||Verify_largest(union_X)){
            ResultSet.push_back(R);
            change_mem(sizeof(int)*R.capacity(), 1);
            return ;
        }
    }
    int left=R.size()%attr_size;
    int wanted=attr_size-left;
    int empty=0;
    for(int i=attr_size-1; i>=attr_size-wanted; i--){
        if(candidates[i].size()==0) empty++;
    }
    if(empty!=0) return;

    vector<int> SwapC;
    for(int i=0; i<R.size(); i++) rvis[R[i]]=1;
    for(int i=0; i<attr_size; i++){
        SwapC.clear();
        for(int j=0; j<candidates[i].size(); j++){
            if(!rvis[candidates[i][j]]){
                SwapC.push_back(candidates[i][j]);
            }
        }
        candidates[i]=SwapC;
    }
    change_mem(sizeof(int)*SwapC.capacity()+sizeof(SwapC), 1);
    #ifdef DEBUG
    printf("C:  ");
    for(int i=0; i<C.size(); i++){
        printf(" %d");
    }
    printf("\n");
    #endif
    for(int i=0; i<R.size(); i++) rvis[R[i]]=0;
    vector<int> newR; newR.clear();
    for(int i=0; i<R.size(); i++)
        newR.push_back(R[i]);
    change_mem(sizeof(newR)+ sizeof(int)*newR.capacity(),1);
    for(int i=0; i<candidates[att_node].size(); i++){
        int cur=candidates[att_node][i];
        newR.push_back(cur);

        vector<int> newC[attr_size]; 
        for(int j=offset[cur]; j<pend[cur]; j++){
            nvis[edge_list[j]]=1;
        }
        for(int j=0; j<attr_size; j++){
            newC[j].clear();
            if(j==att_node){
                for(int k=i+1; k<candidates[j].size(); k++){
                    if(nvis[candidates[j][k]]){
                        newC[j].push_back(candidates[j][k]);
                    }
                }
            }
            else{
                for(int k=0; k<candidates[j].size(); k++){
                    if(nvis[candidates[j][k]]){
                        newC[j].push_back(candidates[j][k]);
                    }
                }
            }
        }
        long long sum_newC=0;
        for(int j=0; j<attr_size; j++) sum_newC+=newC[j].capacity();
        change_mem(sizeof(int)*sum_newC, 1);

        int minval=n;
        for(int j=0; j<attr_size; j++){
            if(newC[j].size()<minval) minval=newC[j].size();
        }
        minval++;
        if(minval*attr_size+newR.size()<threshold*attr_size){
            for(int j=offset[cur]; j<pend[cur]; j++){
                nvis[edge_list[j]]=0;
            }
            newR.pop_back();
            change_mem(sizeof(int)*sum_newC, 0);
            continue;
        }
        /*
        vector<int> newX; newX.clear();
        for(int j=0; j<X.size(); j++){
            if(nvis[X[j]]) newX.push_back(X[j]);
        }
        */
        for(int j=offset[cur]; j<pend[cur]; j++){
            nvis[edge_list[j]]=0;
        }
        vector<int> newX;
        Backtrack_Strong(newR, newC, newX, (att_node+1)%attr_size);
        //X.push_back(cur);
        newR.pop_back();
        change_mem(sizeof(int)*sum_newC, 0);
    }
    change_mem(sizeof(int)*SwapC.capacity()+sizeof(SwapC), 0);
    change_mem(sizeof(newR)+ sizeof(int)*newR.capacity(),0);
}



void Graph::Strong_reduction_Heuristic(){
    colorfulcore();
    printf("left size1=%d\n",left.size());
    int** att_cnt= new int*[max_color];
    for(int i=0; i<max_color; i++) att_cnt[i]= new int[attr_size];
    int* cnt=new int [attr_size];
    int* head= new int[attr_size+1];
    int * nxt=new int[max_color];
    int* vis=new int[n];
    int* fairness_d= new int[n];
    for(int i=0; i<n; i++) vis[i]=1;
    for(int i=0; i<left.size(); i++) vis[left[i]]=0;

    vector<int> upgrade_list;upgrade_list.clear();
    vector<int> upgrade_swap;
    for(int i=0; i<left.size(); i++) upgrade_list.push_back(left[i]);
    while(!upgrade_list.empty()){
        upgrade_swap.clear();
        for(int i=0; i<upgrade_list.size(); i++){
            int nod=upgrade_list[i];
            if(vis[nod]==1) continue;//already deleted
            for(int j=0; j<max_color; j++){
                for(int k=0; k<attr_size; k++) att_cnt[j][k]=0;
            }
            for(int j=offset[nod]; j<pend[nod]; j++){
                int nei=edge_list[j];
                if(vis[nei]==0)
                    att_cnt[color[nei]][attribute[nei]]=1;
            }
            for(int j=0; j<attr_size; j++) cnt[j]=0;//number of colors an attribute owns
            for(int j=0; j<=attr_size; j++) head[j]=n;
            for(int j=0; j<max_color; j++){
                int siz=0, num_;
                for(int k=0; k<attr_size; k++){
                    if(att_cnt[j][k]) {
                        siz++;
                        num_=k;
                    }
                }
                if(siz==1){
                    cnt[num_]++;
                }
                else{
                    nxt[j]=head[siz];
                    head[siz]=j;
                }
            }
            for(int j=1; j<=attr_size; j++){
                for(int color_now=head[j]; color_now!=n; color_now=nxt[color_now]){
                    int minatt=n, minnod;
                    for(int k=0; k<attr_size; k++){
                        if(att_cnt[color_now][k]){
                            if(cnt[k]<minatt){//have the minimum number of color
                                minatt=cnt[k];minnod=k;
                            }
                        }
                    }
                    cnt[minnod]++;
                }
            }
            
            int mincnt=n;
            for(int j=0; j<attr_size; j++){
                if(mincnt>cnt[j]) mincnt=cnt[j];
            }
            fairness_d[nod]=mincnt;
            if(fairness_d[nod]<threshold-1){
                upgrade_swap.push_back(nod);
                vis[nod]=1;
            }
        }
        upgrade_list=upgrade_swap;
    }

    //construct subgraph
    construct_sub_graph(vis);
    printf("left size2=%d\n",left.size());
    delete[] vis;
    printf("reduction done\n");
}

void Graph::Strong_reduction_basic(){
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

    /*recolor
    int* head=new int[n];
    int* nxt= new int[n];
    int* degree=new int[n];
    int* cvis=new int[n];
    int max_degree=0;
    for(int i=0; i<n; i++) degree[i]=0;
    for(int i=0; i<n; i++) head[i]=n;
    for(int i=0; i<n; i++) cvis[i]=0;
    max_color=0;
    for(int i=0; i<left.size(); i++){
        degree[left[i]]=pend[left[i]]-offset[left[i]];
    }
    for(int i=0; i<left.size(); i++){
        int nod=left[i];
        nxt[nod]=head[degree[nod]];
        head[degree[nod]]=nod;
        if(degree[nod]>max_degree)
            max_degree=degree[nod];
    }
    for(int i=0; i<n; i++) color[i]=n;
	max_color = 0;
	for(int ii=max_degree; ii>=1; ii--){
		for(int jj=head[ii]; jj!=n; jj=nxt[jj]){
			int u = jj;
			for(int j = offset[u];j < pend[u];j ++) {
                int c = color[edge_list[j]];
                if(c != n) {
                    cvis[c] = 1;
                }
			}

			for(int j = 0;;j ++){
                if(!cvis[j]) {
                    color[u] = j;
                    if(j > max_color) max_color = j;
                    break;
			    }
            }
			for(int j = offset[u];j < pend[u];j ++) {
                int c = color[edge_list[j]];
                if(c != n) cvis[c] = 0;
			}
		}
	}
    max_color++;
    printf("max_color2=%d\n", max_color);
    */

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
        if(c1<c2){
            if(cboth>(c2-c1)) fairness_d[i]=(cboth-(c2-c1))/2+c2;
            else fairness_d[i]=cboth+c1;
        }
        else if(c1==c2) fairness_d[i]=c1+cboth/2;
        else if(c1>c2){
            if(cboth>(c1-c2)) fairness_d[i]=(cboth-(c1-c2))/2+c1;
            else fairness_d[i]=cboth+c2;
        }
        //printf("fairnessdegree[%d]=%d\n", i, fairness_d[i]);
        //
    }

    //queue<int> Q;
    printf("threshold=%d\n", threshold);
    int deletenode=0;
    int updatenode=0;
    for(int i=0; i<n; i++){
        if(!vis[i]&&fairness_d[i]<threshold){
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
            if(c1<c2){
                if(cboth>(c2-c1)) fairness_d[i]=(cboth-(c2-c1))/2+c2;
                else fairness_d[i]=cboth+c1;
            }
            else if(c1==c2) fairness_d[i]=c1+cboth/2;
            else if(c1>c2){
                if(cboth>(c1-c2)) fairness_d[i]=(cboth-(c1-c2))/2+c1;
                else fairness_d[i]=cboth+c2;
            }
            //printf("fairness degree[%d]=%d\n", i, fairness_d[i]);
            if(fairness_d[i]<threshold){
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
    printf("max mem for strong reduction=%lld\n", max_mem);
    printf("time for second reduction: %s (microseconds)\n", Utility::integer_to_string(tt.elapsed()).c_str());
}

void Graph::Strong_reduction(){
    colorfulcore();  //k->k/2
    int* vis = new int[n];
    for(int i=0; i<n; i++) vis[i]=1;
    for(int i=0; i<left.size(); i++) vis[left[i]]=0;

    int* fairness_d=new int[n];
    for(int i=0; i<n; i++) fairness_d[i]=0;
    //calculate groups
    int** cntGroup=new int*[n];
    for(int i=0; i<n; i++){
        cntGroup[i]= new int[3];
        for(int j=0; j<3; j++){
            cntGroup[i][j]=0;
        }
    }


    int** colorful_d=new int*[n];
    int** c_record= new int*[n];
    for(int i=0; i<n; i++){
        c_record[i] = new int[3];
        for(int j=0; j<3; j++){
            c_record[i][j]=0;
        }
    }
    for(int i=0; i<n; i++){
        colorful_d[i]= new int[max_color];
        for(int j=0; j<max_color; j++){
            colorful_d[i][j]=0;
        }
    }

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
                    colorful_d[i][j]=num*10+0;
                }
                else{
                    c2++;
                    colorful_d[i][j]=num*10+1;
                }
            }
            else if(ccnt>1){
                cboth++;
                colorful_d[i][j]=2;
            }
        }

        //calculate fairness degree
        c_record[i][1]=c1;
        c_record[i][2]=c2;
        c_record[i][3]=cboth;
        if(c1<c2){
            if(cboth>=(c2-c1)) fairness_d[i]=(cboth-(c2-c1))/2+c2;
            else fairness_d[i]=cboth+c1;
        }
        else if(c1==c2) fairness_d[i]=c1+c2+cboth/2;
        else if(c1>c2){
            if(cboth>=(c1-c2)) fairness_d[i]=(cboth-(c1-c2))/2+c2;
            else fairness_d[i]=cboth+c2;
        }
        //
    }

    //reduce graph  improved type
    
    queue<int> Q;
    int deletenode=0;
    for(int i=0; i<n; i++){
        if(!vis[i]&&fairness_d[i]<threshold-1){
            Q.push(i);
            vis[i]=1;
        }
    }
    int num;
    while(!Q.empty()){
        int cur=Q.front();
        Q.pop();
        for(int i=offset[cur]; i<pend[cur]; i++){
            int neighbor=edge_list[i];
            if(!vis[neighbor]){
                typ=colorful_d[neighbor][color[cur]]%10;
                if(typ<2){
                    num=colorful_d[neighbor][color[cur]]/10;
                    num--;
                    if(num){
                        colorful_d[neighbor][color[cur]]=num*10+typ;
                    }
                    else{
                        c_record[neighbor][typ]--;
                    }
                }
                else{//if type==cboth, check if the type is changed
                    typ=-1;num=0;
                    int flag=1;
                    for(int j=offset[neighbor]; j<pend[neighbor]; j++){
                        if(edge_list[j]!=cur&&color[edge_list[j]]==color[cur]){
                            num++;
                            if(typ!=-1&&attribute[edge_list[j]]!=typ){
                                flag=0;break;
                            }
                            typ=attribute[edge_list[j]];
                        }
                    }
                    if(flag==1){//become c1 or c2
                        colorful_d[neighbor][color[cur]]=num*10+typ;
                        c_record[neighbor][typ]++;
                        c_record[neighbor][2]--;
                    }
                }
                //re-calculate fairness degree
                if(fairness_d[neighbor]<threshold-1){
                    Q.push(i);
                    vis[neighbor]=1;
                }
            }
        }
    }
    left.clear();
    for(int i=0; i<n; i++){
        if(!vis[i])
            left.push_back(i);
    }
    construct_sub_graph(vis);
    delete[] c_record;
    delete[] colorful_d;
}

//fairness degree
//不管颜色，只看数量
vector<int> Graph::GetFairnessOrdering(int att_siz){
    //printf("get fairness ordering\n");
    int** d;//fairness degree
    change_mem(sizeof(int)*n*attr_size,1);
    int* min_d;//minimum among all attributes
    change_mem(sizeof(int)*n, 1);
    d=new int*[n];
    min_d=new int[n];
    int* vis=new int[n];
    change_mem(sizeof(int)*n, 1);
    for(int i=0; i<n; i++) min_d[i]=n;
    for(int i=0; i<n; i++){
        d[i]=new int[att_siz];
        for(int j=0; j<att_siz; j++)
            d[i][j]=0;
    }
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<component.size(); i++) vis[component[i]]=1;
    
    for(int i=0; i<component.size(); i++){
        for(int j=offset[component[i]]; j<pend[component[i]]; j++){
            if(vis[edge_list[j]]){
               // printf("yes\n");
                int a_nei=attribute[edge_list[j]];//attribute of neighbor
                d[component[i]][a_nei]++;
            }
        }
    }
    
    int min_degree_node, min_degree=n;
    for(int i=0; i<n; i++){
        if(!vis[i]) continue;
        for(int j=0; j<att_siz; j++){
            if(d[i][j]<min_d[i]) min_d[i]=d[i][j];
        }
    }
    
    int* component_idx=new int[n];
    change_mem(sizeof(int)*n, 1);
    for(int i=0; i<component.size(); i++)
        component_idx[component[i]]=i;
    vector<int> peeling_order; peeling_order.clear();
    int* degree_arr=new int[component.size()];
    change_mem(sizeof(int)*n, 1);
    for(int i=0; i<component.size(); i++)
        degree_arr[i]=min_d[component[i]];
    /*
    for(int i=0; i<component.size(); i++){
        printf("degree_arr[%d]=%d\n", i, degree_arr[i]);
    }
    */
    int* ordered_c=new int[component.size()];
    change_mem(sizeof(int)*component.size(), 1);
    for(int i=0; i<component.size();i++) ordered_c[i]=i;
    ListLinearHeap *heap = new ListLinearHeap(component.size(), component.size());
    change_mem(heap->get_max_occupy(),1);
    heap->init(component.size(), component.size(), ordered_c, degree_arr);
    vector<int> reversed_peeling_order;
    for(int i=0; i<component.size(); i++){
        int nod, key;
        heap->pop_min(nod, key);
        nod=component[nod];
        reversed_peeling_order.push_back(nod);
        vis[nod]=0;
        for(int j=offset[nod]; j<pend[nod]; j++){
            if(vis[edge_list[j]]==1){
                if((--d[edge_list[j]][attribute[nod]])<min_d[edge_list[j]]){
                    int off=min_d[edge_list[j]]-d[edge_list[j]][attribute[nod]];
                    min_d[edge_list[j]]=d[edge_list[j]][attribute[nod]];
                    heap->decrement(component_idx[edge_list[j]], off);
                }
            }
        }
    }
    for(int i=0; i<component.size(); i++)
        peeling_order.push_back(reversed_peeling_order[i]);
    change_mem(sizeof(peeling_order)+sizeof(int)*peeling_order.capacity(), 1);
    change_mem(sizeof(reversed_peeling_order)+sizeof(int)*reversed_peeling_order.capacity(), 1);
    
    
    change_mem(sizeof(peeling_order)+sizeof(int)*peeling_order.capacity(), 0);
    change_mem(sizeof(reversed_peeling_order)+sizeof(int)*reversed_peeling_order.capacity(), 0);
    change_mem(heap->get_max_occupy(),0);  
    delete[] d;     change_mem(sizeof(int)*n*attr_size, 0);
    delete[] min_d; change_mem(sizeof(int)*n, 0);
    delete[] vis;   change_mem(sizeof(int)*n, 0);
    delete[] component_idx;   change_mem(sizeof(int)*n, 0);
    delete[] degree_arr;      change_mem(sizeof(int)*n, 0);
    delete[] ordered_c;       change_mem(sizeof(int)*n, 0);
    //printf("compute done!\n");
    return peeling_order;
}

vector<int> Graph::GetColorfulFairnessOrdering(){
    //printf("get colorful fairness ordering\n");
    clear_local_mem();
    change_local_mem(sizeof(int)*n,1);
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
    change_local_mem(sizeof(int)*n*3,1);
    int* queue_n= new int[component.size()];
    change_local_mem(sizeof(int)*component.size(),1);
    for(int i=0; i<component.size(); i++) queue_n[i]=i;
    int* queue_fair= new int[component.size()];
    change_local_mem(sizeof(int)*component.size(),1);
    for(int i=0; i<component.size(); i++) queue_fair[i]=fairness_d[component[i]];
    int* idx= new int[n];
    change_local_mem(sizeof(int)*n,1);
    for(int i=0; i<component.size(); i++)
        idx[component[i]]=i;
       printf("step3");

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
            if(c1<c2){
                if(cboth>(c2-c1)) fairness_d[nei]=(cboth-(c2-c1))/2+c2;
                else fairness_d[nei]=cboth+c1;
            }
            else if(c1==c2) fairness_d[nei]=c1+cboth/2;
            else if(c1>c2){
                if(cboth>(c1-c2)) fairness_d[nei]=(cboth-(c1-c2))/2+c1;
                else fairness_d[nei]=cboth+c2;
            }

            heap->decrement(idx[nei], old_fairness-fairness_d[nei]);//decrease the degree by 1
        }
    }
    change_local_mem(sizeof(peeling_order)+sizeof(int)*peeling_order.capacity(), 1);
    delete heap;
    change_local_mem(heap->get_max_occupy(), 0);
    delete[] cntGroup;
    change_local_mem(sizeof(int)*n*3,0);
    delete[] vis;
    change_local_mem(sizeof(int)*n, 0);
    delete[] queue_n;
    change_local_mem(sizeof(int)*component.size(),0);
    delete[] queue_fair;
    change_local_mem(sizeof(int)*component.size(),0);
    delete[] idx;
    change_local_mem(sizeof(int)*n,0);
    change_local_mem(sizeof(peeling_order)+sizeof(int)*peeling_order.capacity(), 0);
    //printf("compute done!\n");
    return peeling_order;
}

vector<int> Graph::GetColorfulFairnessOrdering_Heuristic(){
    set<int> c_attr[max_color];
    set<int>::iterator it;
    set<pair<int, int> > color_idx[max_color];
    int* colorful_degree_heu= new int[n];
    for(int i=0; i<n; i++) colorful_degree_heu[i]=0;
    int* attr_count = new int[attr_size];
    set<pair<int, int> > ordered_count;
    set<pair<int, int> >:: iterator cnt_it;
    vector<int> cnt[attr_size];

    for(int i=0; i<component.size(); i++){
        for(int j=0; j<max_color; j++) c_attr[j].clear();
        for(int j=offset[component[i]]; j<pend[component[i]]; j++){
            c_attr[color[edge_list[j]]].insert(attribute[edge_list[j]]);
        }
        for(int j=0; j<attr_size; j++) attr_count[j]=0;
        for(int j=0; j<max_color; j++) color_idx[j].clear();
        for(int j=0; j<max_color; j++){
            for(it=c_attr[j].begin(); it!=c_attr[j].end(); it++)
                color_idx[*it].insert(make_pair(c_attr[j].size(), j));
        }
        for(int j=0; j<attr_size; j++) cnt[j].clear();
        for(int j=0; j<attr_size; j++){
            for(cnt_it=color_idx[j].begin(); cnt_it!=color_idx[j].end(); cnt_it++){
                cnt[j].push_back(cnt_it->first);
            }
        }
        int min_sum=n;
        for(int j=0;j<attr_size; j++){
            double sum=0;
            for(int k=0; k<cnt[j].size(); k++){
                sum+=1.0/cnt[j][k];
            }
            if((int)sum <min_sum) min_sum=sum;
        }
        colorful_degree_heu[component[i]]=min_sum;
    }
    vector<int> peeling_order;peeling_order.clear();
    int* ordered_c=new int[component.size()];
    for(int i=0; i<component.size();i++) ordered_c[i]=i;
    int* degree_arr=new int[component.size()];
    for(int i=0; i<component.size(); i++) degree_arr[i]=colorful_degree_heu[component[i]];
    
    ListLinearHeap *heap = new ListLinearHeap(component.size(), component.size());
    heap->init(component.size(), component.size(), ordered_c, degree_arr);
    vector<int> reversed_peeling_order;
    for(int i=0; i<component.size(); i++){
        int nod, key;
        heap->pop_min(nod, key);
        nod=component[nod];
        reversed_peeling_order.push_back(nod);
    }
    for(int i=0; i<component.size(); i++)
        peeling_order.push_back(reversed_peeling_order[i]);
    return peeling_order;
}

vector<int> Graph::GetColorfulFairnessOrdering_Heuristic_2(){
 
    clear_local_mem();
    int* head;
    int* nxt;

    int** att_cnt= new int*[max_color];
    for(int i=0; i<max_color; i++) att_cnt[i]= new int[attr_size];
    int* cnt=new int [attr_size];
    head= new int[attr_size+1];
    nxt=new int[max_color];
    change_local_mem(sizeof(int)*max_color*attr_size, 1);
    change_local_mem(sizeof(int)*attr_size, 1);
    change_local_mem(sizeof(int)*(attr_size+1), 1);
    change_local_mem(sizeof(int)*max_color, 1);

    int* fairness_d= new int[n];
    change_local_mem(sizeof(int)*n, 1);
    int* vis=new int[n];
    change_local_mem(sizeof(int)*n, 1);
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<component.size(); i++) vis[component[i]]=1;

    vector<int> peeling_order;
    int nodecnt=0;
    vector<int> upgrade_list;upgrade_list.clear();
    for(int i=0; i<component.size(); i++) upgrade_list.push_back(component[i]);
    while(nodecnt<component.size()){
        for(int i=0; i<upgrade_list.size(); i++){
            int nod=upgrade_list[i];
            if(vis[nod]==0) continue;
            for(int j=0; j<max_color; j++){
                for(int k=0; k<attr_size; k++) att_cnt[j][k]=0;
            }
            for(int j=offset[nod]; j<pend[nod]; j++){
                int nei=edge_list[j];
                if(vis[nei]==1)
                    att_cnt[color[nei]][attribute[nei]]=1;
            }
            for(int j=0; j<attr_size; j++) cnt[j]=0;//number of colors an attribute owns
            for(int j=0; j<=attr_size; j++) head[j]=n;
            for(int j=0; j<max_color; j++){
                int siz=0, num_;
                for(int k=0; k<attr_size; k++){
                    if(att_cnt[j][k]) {
                        siz++;
                        num_=k;
                    }
                }
                if(siz==1){
                    cnt[num_]++;
                }
                else{
                    nxt[j]=head[siz];
                    head[siz]=j;
                }
            }
           
            for(int j=1; j<=attr_size; j++){
                for(int color_now=head[j]; color_now!=n; color_now=nxt[color_now]){
                    int minatt=n, minnod;
                    for(int k=0; k<attr_size; k++){
                        if(att_cnt[color_now][k]){
                            if(cnt[k]<minatt){//have the minimum number of color
                                minatt=cnt[k];minnod=k;
                            }
                        }
                    }
                    cnt[minnod]++;
                }
            }
            
            int mincnt=n;
            for(int j=0; j<attr_size; j++){
                if(mincnt>cnt[j]) mincnt=cnt[j];
            }
            fairness_d[nod]=mincnt;
        }
        int min_fairness=n, min_nod;
        for(int j=0; j<component.size(); j++){
            if(vis[component[j]]&&min_fairness>fairness_d[component[j]]){
                min_fairness=fairness_d[component[j]];
                min_nod=component[j];
            }
        }
        peeling_order.push_back(min_nod);
        vis[min_nod]=0;
        nodecnt++;
        upgrade_list.clear();
        for(int j=offset[min_nod]; j<pend[min_nod]; j++){
            if(vis[edge_list[j]]) upgrade_list.push_back(edge_list[j]);
        }
    }
    change_local_mem(sizeof(upgrade_list)+sizeof(int)*upgrade_list.capacity(), 1);
    change_local_mem(sizeof(peeling_order)+sizeof(int)*peeling_order.capacity(), 1);

    change_local_mem(sizeof(upgrade_list)+sizeof(int)*upgrade_list.capacity(), 0);
    change_local_mem(sizeof(peeling_order)+sizeof(int)*peeling_order.capacity(), 0);
    delete[] head;change_local_mem(sizeof(int)*(attr_size+1), 0);
    delete[] nxt; change_local_mem(sizeof(int)*max_color, 0);
    delete[] att_cnt;change_local_mem(sizeof(int)*max_color*attr_size, 0);
    delete[] cnt;change_local_mem(sizeof(int)*attr_size, 0);
    delete[] vis;change_local_mem(sizeof(int)*n, 0);

    return peeling_order;
}

vector<int> Graph::core_decomposition(int target_attr){
    printf("target attr=%d\n", target_attr);
    int** core_level=new int*[n];
    for(int i =0; i<n; i++){
        core_level[i]=new int[attr_size+1];
        for(int j=0; j<=attr_size; j++)
            core_level[i][j]=0;
    }
    vector<int> comp_special;
    int lop=0;
    int reduce_cnt=0;
    vector<int> delete_set;delete_set.clear();
    while(lop<attr_size){
        comp_special.clear();
        for(int i=0; i<component.size(); i++){
            if(target_attr+lop<attr_size){
                if(attribute[component[i]]>=target_attr&&
                attribute[component[i]]<=target_attr+lop) 
                    comp_special.push_back(component[i]);
            }
            else{
                if(attribute[component[i]]>=target_attr&&attribute[component[i]]<=target_attr+lop
                || attribute[component[i]]>=0&&attribute[component[i]]<=(target_attr+lop)%attr_size)
                     comp_special.push_back(component[i]);
            }
        }
        int* degree = new int[comp_special.size()];
        int* queue_n= new int[comp_special.size()];
        int* idx= new int[n];
        int* vis= new int[n];
        for(int i=0; i<n; i++) vis[i]=0;
        for(int i=0; i<comp_special.size(); i++) degree[i]=0;
        for(int i=0; i<comp_special.size(); i++) vis[comp_special[i]]=1;
        for(int i=0; i<comp_special.size(); i++){
            for(int j=offset[comp_special[i]]; j<pend[comp_special[i]]; j++)
                if(vis[edge_list[j]]) degree[i]++;
            queue_n[i]=i;
            idx[comp_special[i]]=i;
        }
       // printf("computing...com.size=%d\n", comp_special.size());
        ListLinearHeap *heap = new ListLinearHeap(comp_special.size(), comp_special.size());//nodes, max degree
        heap->init(comp_special.size(), comp_special.size(), queue_n, degree);//number of nodes, max degree, start position of arr, start pos of degreesa
        
        vector<ui> res;
        int max_core=0;
        for(ui i = 0;i < comp_special.size();i ++) {
            ui u, key;
            heap->pop_min(u, key);
            int node=comp_special[u];
            if(key>max_core) max_core=key;
            core_level[node][lop+1]=max_core;
            //ordered_arr.push_back(node);
            vis[node] = 0;
            for(ui j = offset[node];j < pend[node];j ++) if(vis[edge_list[j]] == 1) {
                heap->decrement(idx[edge_list[j]], 1);//decrease the degree by 1
            }
        }
        delete heap;
        delete[] degree;
        delete[] queue_n;
        delete[] idx;
        delete[] vis;
        lop++;
    }
    int* min_core= new int[n];
    for(int i=0; i<n; i++) min_core[i]=n;
    for(int i=0; i<n ;i++) core_level[i][0]=-1;
    vector<int> component_reduce; component_reduce.clear();
    for(int i=0; i<component.size(); i++){
        int nod=component[i];
        if(attribute[nod]!=target_attr) continue;
        component_reduce.push_back(nod);
        for(int j=1; j<=attr_size; j++){
            if(min_core[nod]>core_level[nod][j]-core_level[nod][j-1]){
                min_core[nod]=core_level[nod][j]-core_level[nod][j-1];
            }
        }  
    }
    int* queue_n= new int[component_reduce.size()];
    int* degree= new int[component_reduce.size()];
    for(int i=0; i<component_reduce.size(); i++) queue_n[i]=i;
    for(int i=0; i<component_reduce.size(); i++) degree[i]=min_core[component_reduce[i]];
   // for(int i=0; i<component.size(); i++) printf("degree[%d]=%d\n", i, degree[i]);
    ListLinearHeap *heap = new ListLinearHeap(component_reduce.size(), component_reduce.size());//nodes, max degree
    heap->init(component_reduce.size(), component_reduce.size(), queue_n, degree);//number of nodes, max degree, start position of arr, start pos of degreesa
    vector<int> ordered;ordered.clear();
    for(ui i = 0;i < component_reduce.size();i ++) {
        ui u, key;
        heap->pop_min(u, key);
        ordered.push_back(component_reduce[u]);
    }
    delete heap;
    return ordered;
}

vector<int> Graph::core_decomposition_cut(){
    vector<int> comp_special;
    int lop=0;
    int reduce_cnt=0;
    vector<int> delete_set;delete_set.clear();
    while(lop<attr_size){
        comp_special.clear();
        for(int i=0; i<left.size(); i++){
            if(attribute[left[i]]==lop) comp_special.push_back(left[i]);
        }
        int* degree = new int[comp_special.size()];
        int* queue_n= new int[comp_special.size()];
        int* idx= new int[n];
        int* vis= new int[n];
        for(int i=0; i<n; i++) vis[i]=0;
        for(int i=0; i<comp_special.size(); i++) degree[i]=0;
        for(int i=0; i<comp_special.size(); i++) vis[comp_special[i]]=1;
        for(int i=0; i<comp_special.size(); i++){
            for(int j=offset[comp_special[i]]; j<pend[comp_special[i]]; j++)
                if(vis[edge_list[j]]) degree[i]++;
            queue_n[i]=i;
            idx[comp_special[i]]=i;
        }
        printf("computing...com.size=%d\n", comp_special.size());
        ListLinearHeap *heap = new ListLinearHeap(comp_special.size(), comp_special.size());//nodes, max degree
        heap->init(comp_special.size(), comp_special.size(), queue_n, degree);//number of nodes, max degree, start position of arr, start pos of degreesa
        
        vector<ui> res;
        int max_core=0;
        for(ui i = 0;i < comp_special.size();i ++) {
            ui u, key;
            heap->pop_min(u, key);
            int node=comp_special[u];
            if(key>max_core) max_core=key;
            if(max_color<threshold-1) {
                reduce_cnt++;
                delete_set.push_back(node);
            }
            //ordered_arr.push_back(node);
            vis[node] = 0;
            for(ui j = offset[node];j < pend[node];j ++) if(vis[edge_list[j]] == 1) {
                heap->decrement(idx[edge_list[j]], 1);//decrease the degree by 1
            }
        }
        
        delete heap;
        delete[] degree;
        delete[] queue_n;
        delete[] idx;
        delete[] vis;
        lop++;
    }
    printf("reduce_cnt=%d\n", reduce_cnt);
    return delete_set;
}

long long Graph::get_max(){
    return this->max_mem;
}

vector<int> Graph::GetColorfulOrdering(){
    int* vis=new int[n];
    for(int i=0; i<n; i++) vis[i]=0;
    for(int i=0; i<component.size(); i++) vis[component[i]]=1;

    int*** d;//fairness degree
    int** r;

    int* min_d;//minimum among all attributes
    int* component_idx=new int[n];
    for(int i=0; i<component.size(); i++)
        component_idx[component[i]]=i;
    min_d=new int[n];
    for(int i=0; i<n; i++) min_d[i]=n;
    r=new int*[component.size()];
    for(int i=0; i<component.size(); i++){
        r[i]= new int[attr_size];
        for(int j=0; j<attr_size; j++)
            r[i][j]=0;
    }
    d=new int** [component.size()];
    for(int i=0; i<component.size(); i++){
        d[i] =new int*[attr_size];
        for(int j=0; j<attr_size; j++){
            d[i][j]= new int [max_color];
            for(int k=0; k<max_color; k++)
                d[i][j][k]=0;
        }
    }

    for(int i=0; i<component.size(); i++){
        int cnow=component[i];
        for(int j=offset[cnow]; j<pend[cnow]; j++){
            int nei=edge_list[j];
            if(vis[nei]==1)
                if((d[i][attribute[nei]][color[nei]]++) ==0)
                    r[i][attribute[nei]]++;
        }

        for(int j=0; j<attr_size; j++){
            if(r[i][j]<min_d[i]) min_d[i]=r[i][j];
        }
    }

    int* degree_arr=new int[component.size()];

    for(int i=0; i<component.size(); i++)
        degree_arr[i]=min_d[i];
    int* ordered_c=new int[component.size()];
    for(int i=0; i<component.size();i++) ordered_c[i]=i;
    idx_pos = new int[n];
    int cnt=0;
    ListLinearHeap *heap = new ListLinearHeap(component.size(), component.size());
    heap->init(component.size(), component.size(), ordered_c, degree_arr);
    vector<int> reversed_peeling_order;reversed_peeling_order.clear();
    for(int i=0; i<component.size(); i++){
        int nod, key;
        heap->pop_min(nod, key);
        nod=component[nod];
        //printf("nod=%d, key=%d\n", nod, key);
        reversed_peeling_order.push_back(nod);
        idx_pos[nod]= cnt ++;
        vis[nod]=0;
        for(int j=offset[nod]; j<pend[nod]; j++){
            if(vis[edge_list[j]]==1){
                int idx=component_idx[edge_list[j]];
                if(--d[idx][attribute[nod]][color[nod]] <=0){
                    r[idx][attribute[nod]]--;
                    if(min_d[idx]>r[idx][attribute[nod]]){
                        int off =min_d[idx]- r[idx][attribute[nod]];
                        min_d[idx]=r[idx][attribute[nod]];
                        heap->decrement(idx, off);
                    }
                }
            }
        }
    }



    delete[] d;
    delete[] r;
    delete[] min_d;
    delete[] vis;
    delete[] component_idx;
    delete[] degree_arr;
    delete[] ordered_c;    

    return reversed_peeling_order;
}

vector<int> Graph::GetDegreeOrdering(){
    int* d= new int[n];
    int* head=new int [n];
    int* nxt= new int[n];
    for(int i=0; i<n; i++) head[i]=n;
    int max_d=0;
    for(int i=0; i<component.size(); i++){
        int cur=component[i];
        d[cur]=pend[cur]- offset[cur];
        nxt[cur]=head[d[cur]];
        head[d[cur]]=cur;
        if(d[cur]>max_d) max_d=d[cur];
    }
    vector<int> ordering;ordering.clear();
    for(int i=0; i<=max_d; i++){
        for(int j=head[i]; j!=n; j=nxt[j]){
            ordering.push_back(j);
        }
    }
    delete[] d;
    delete[] head;
    delete[] nxt;
    return ordering;
}

void Graph::Strong_reduction_Three(){
    colorfulcore();
    Timer tt;

    int* vis = new int[n];
    for(int i=0; i<n; i++) vis[i]=1;
    for(int i=0; i<left.size(); i++) vis[left[i]]=0;

    fairness_d=new int[n];
    for(int i=0; i<n; i++) fairness_d[i]=0;
    //calculate groups
    int** cntGroup=new int*[n];
    for(int i=0; i<n; i++){
        cntGroup[i]= new int[attr_size];
        for(int j=0; j<attr_size; j++){
            cntGroup[i][j]=0;
        }
    }

    int cboth, ccnt, typ , call, sum_res;
    int* c_list= new int[attr_size+1];
    int* c_assist=new int[attr_size+1];
    int totcnt=0;
    for(int i=0; i<n; i++){
        if(vis[i]) continue;
        cboth=0;call=0;sum_res=0;
        for(int j=0; j<attr_size+1; j++)
            c_list[j]=c_assist[j]=0;
        for(int j=0; j<=max_color; j++){
            for(int k=0; k<attr_size; k++)
                cntGroup[j][k]=0;
        }
        for(int j=offset[i]; j<pend[i]; j++){
            if(!vis[edge_list[j]])
                ++cntGroup[color[edge_list[j]]][attribute[edge_list[j]]];
        }
        for(int j=0; j<max_color; j++){
            ccnt=0;typ=0;sum_res=0;
            for(int k=attr_size-1; k>=0; k--){
                if(cntGroup[j][k]){
                    ccnt++;
                    typ=k;
                    sum_res+=k;
                }
            }
            if(ccnt==1)
                c_list[typ]++;
            else if(ccnt<=2){
                if(sum_res==1)
                    cboth++;
                else {
                    c_list[typ]++;
                    c_assist[typ]++;
                }
            }
            else if(ccnt>2){
                call++;
            }
        }
        c_list[attribute[i]]++;
        if(c_list[0]<c_list[1]){
            if(cboth>(c_list[1]-c_list[0])){
                int tmp=cboth-(c_list[1]-c_list[0]);
                c_list[1]=tmp/2+c_list[1];
                c_list[0]=c_list[1];
                if(tmp%2==1) c_list[0]++;
            }
            else{
                c_list[0]+=cboth;
            }
        }
        else if(c_list[0]==c_list[1]) {
            c_list[0]+=cboth/2;
            c_list[1]+=cboth/2;
            if(cboth%2==1) c_list[0]++;
        }
        else if(c_list[0]>c_list[1]){
            if(cboth>(c_list[0]-c_list[1])){
                int tmp=cboth-(c_list[0]-c_list[1]);
                c_list[0]=tmp/2+c_list[0];
                c_list[1]=c_list[0];
                if(tmp%2==1) c_list[0]++;
            } 
            else{
                c_list[1]+=cboth;
            }
        }
        int lack=0;

        for(int j=0; j<attr_size-1; j++){
            if(c_list[j]<threshold){
                call-=(threshold-c_list[j]);
                if(call<0){
                    vis[i]=1;totcnt++;
                }
            }
        }
        if(vis[i]==0){
            if(call>=0){
                c_list[attr_size-1]+=call;
            }
            if(c_list[attr_size-1]<threshold){
                for(int j=0; j<attr_size-1; j++){
                    c_list[attr_size-1]+=min(c_list[j]-threshold,c_assist[j]);
                }
                if(c_list[attr_size-1]<threshold){
                    vis[i]=1;totcnt++;
                }
            }
        }
    }
    printf("totcnt=%d\n", totcnt);
    construct_sub_graph(vis);
    //#ifdef DEBUG
    ofstream leftofs("left.txt");
    for(int i=0; i<left.size(); i++){
        leftofs<<left[i]<<"\n";
    }
    //#endif
}
