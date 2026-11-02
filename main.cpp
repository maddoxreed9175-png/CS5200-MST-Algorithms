#include <iostream>


#include "kruskal.h"
#include "prim.h"
#include "generate.h"

using namespace std;

int main(){
    int V = 3; // Vertices
    int E = 2; // Edges

    generateConnectedGraph(V, E);
    printGraph();
    MSTResult res = kruskal(adj);
    printMST(res);
    res = prim(adj);
    printMST(res);

    return 0;
}
