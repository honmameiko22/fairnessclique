#include "Graph.h"
#include "Utility.h"
using namespace std;

int main(int argc, char *argv[]){
    if(argc<5){
        printf("enter 5 paramaters at least!\n");
        return 0;
    }
    printf("Graph=%s attribute=%s algorithm=%s threshold=%s order=%s delta=%s\n", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
    Graph* graph = new Graph(argv[1]);
    graph->ReadGraph(argv[1], argv[2], atoi(argv[4]), "r");//read and store graph
    if(strcmp(argv[3],"FairClique")==0) graph->FindClique(argv[5]);
    if(strcmp(argv[3],"StrongClique")==0 ) graph->FindStrongClique(argv[5]);
    //based on alternative selection
    if(strcmp(argv[3],"RelativeStrong") ==0 ) graph->FindRelatedFairClique_S(argv[5], atoi(argv[6]));
    //base on weak fair clique enumeration
    if(strcmp(argv[3],"RelativeWeak") ==0 ) graph->FindRelatedFairClique_W("colorful", atoi(argv[6]));
    if(strcmp(argv[3],"Baseline")==0 ) graph->Baseline(argv[5],atoi(argv[6])); //order strong weak
    printf("algorithm done\n");
    printf("max_mem=%lld\n", graph->get_max());
    return 0;
}