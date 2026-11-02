#include <iostream>
#include <vector>
#include <random>
#include <set>
#include <algorithm>

using namespace std;

using Graph = vector<vector<pair<int, int>>>;

// Adjacency list: vector of pairs {neighbor, weight}
vector<vector<pair<int, int>>> adj;
int N, M;

// Helper to get random edge weight
int getWeight(mt19937& gen) {
    uniform_int_distribution<int> dist(1, 100);
    return dist(gen);
}

// Helper to ensure unique undirected edges
bool addEdge(int u, int v, int w, mt19937& gen, set<pair<int, int>>& edges) {
    if (u == v) return false; // No self-loops
    pair<int, int> edge = {min(u, v), max(u, v)};
    if (edges.count(edge)) return false; // Already exists

    adj[u].push_back({v, w});
    adj[v].push_back({u, w});
    edges.insert(edge);
    return true;
}

// Generate a connected graph
Graph generateConnectedGraph(int n, int m) {
    N = n;
    M = m;
    adj.assign(N, vector<pair<int, int>>());
    
    mt19937 gen(random_device{}());
    set<pair<int, int>> edges;

    // 1. Build Spanning Tree (Guarantees Connectivity)
    // Simple approach: Connect vertex i to a random previous vertex i-1 -> 0
    // Better: Random Prim's or just connect 0-1, 1-2... then shuffle
    // Here we use a simple chain + shuffle for simplicity, or random spanning tree logic.
    // For a robust random spanning tree, one might use Kruskal's with random weights.
    
    // Let's use a simple connected seed: connect 0-1, 1-2, ..., (N-2)-(N-1)
    for (int i = 0; i < N - 1; i++) {
        int u = i;
        int v = i + 1;
        int w = getWeight(gen);
        addEdge(u, v, w, gen, edges);
    }

    // 2. Add Remaining Edges
    // We need M - (N-1) more edges
    int edgesNeeded = M - (N - 1);
    
    // If M < N-1, it's impossible to be connected
    if (edgesNeeded < 0) {
        cout << "Error: Cannot have fewer edges than N-1 for a connected graph." << endl;
        return adj;
    }

    while (edgesNeeded > 0) {
        int u = uniform_int_distribution<int>(0, N - 1)(gen);
        int v = uniform_int_distribution<int>(0, N - 1)(gen);
        int w = getWeight(gen);
        
        if (addEdge(u, v, w, gen, edges)) {
            edgesNeeded--;
        }
    }

    return adj;
}

void printGraph() {
    for (int i = 0; i < N; i++) {
        cout << i << ": ";
        for (auto& edge : adj[i]) {
            cout << "(" << edge.first << ", " << edge.second << ") ";
        }
        cout << endl;
    }
}
