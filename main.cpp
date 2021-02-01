#include "Graph.h"
#include "Utility.h"
using namespace std;

//argv[1]->input graph  argv[2]->attribute 
//argv[3]->algorithm    argv[4]->threshold;
//argv[5]->order
//(option)argv[6]->verify data
int main(int argc, char *argv[]){
    if(argc<5){
        printf("enter 5 paramaters!\n");
        return 0;
    }
    printf("Graph=%s attribute=%s algorithm=%s threshold=%s order=%s\n", argv[1], argv[2], argv[3], argv[4], argv[5]);
    Graph* graph = new Graph(argv[1]);
    graph->ReadGraph(argv[1], argv[2], atoi(argv[4]), "r");//read and store graph
    if(strcmp(argv[3],"ColorfulDegree")==0) graph->colorfulcore();
    if(strcmp(argv[3],"FairClique")==0) graph->FindClique(argv[5]);
    if(strcmp(argv[3],"StrongClique")==0 ) graph->FindStrongClique(argv[5]);
    if(strcmp(argv[3],"VerifyW")==0) graph->VerifyWClique(argv[6]);
    if(strcmp(argv[3],"VerifyS")==0) graph->VerifySClique(argv[6]);
    printf("algorithm done\n");
    printf("max_mem=%lld\n", graph->get_max());
    return 0;
}