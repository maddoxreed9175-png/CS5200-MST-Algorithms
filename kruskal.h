#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

// Disjoint Set Union with path compression and union by size.
class DisjointSet {
public:
    explicit DisjointSet(int n) : parent(n), sz(n, 1), components(n) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    // Find with path compression (iterative, avoids deep recursion).
    int find(int x) {
        int root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != root) {
            int next = parent[x];
            parent[x] = root;
            x = next;
        }
        return root;
    }

    // Union by size. Returns true if a merge happened, false if already joined.
    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (sz[a] < sz[b]) std::swap(a, b);
        parent[b] = a;
        sz[a] += sz[b];
        --components;
        return true;
    }

    bool connected(int a, int b) { return find(a) == find(b); }
    int componentCount() const { return components; }

private:
    std::vector<int> parent;
    std::vector<int> sz;
    int components;
};

// Shared types (also defined in prim.h; guarded so both headers can be included together).
#ifndef MST_TYPES_DEFINED
#define MST_TYPES_DEFINED

struct Edge {
    int u, v, w;
};

struct MSTResult {
    std::vector<Edge> edges;  // edges chosen for the MST (or spanning forest)
    long long totalWeight = 0;
    bool connected = true;    // false if the graph is disconnected (result is a forest)
};


// Prints the MST edges and total weight to the given stream (default: std::cout).
// Output order is deterministic regardless of algorithm: each edge is printed with
// its smaller vertex first, and edges are sorted by (smaller vertex, larger vertex,
// weight), so the list starts with the lowest-numbered vertex.
inline void printMST(const MSTResult& res, std::ostream& out = std::cout) {
    std::vector<Edge> sorted = res.edges;  // copy; res is left untouched
    for (Edge& e : sorted) {
        if (e.u > e.v) std::swap(e.u, e.v);
    }
    std::sort(sorted.begin(), sorted.end(), [](const Edge& a, const Edge& b) {
        if (a.u != b.u) return a.u < b.u;
        if (a.v != b.v) return a.v < b.v;
        return a.w < b.w;
    });

    out << "MST edges:\n";
    for (const Edge& e : sorted) {
        out << e.u << " - " << e.v << "  (weight " << e.w << ")\n";
    }
    out << "Total weight: " << res.totalWeight << "\n";
    if (!res.connected) {
        out << "Note: graph is disconnected, so this is a spanning forest.\n";
    }
}

#endif // MST_TYPES_DEFINED

// Kruskal's algorithm.
// adj[u] holds pairs (v, weight) for every edge u -- v.
// For undirected graphs each edge normally appears in both adj[u] and adj[v];
// duplicates are harmless since the DSU rejects the second copy.
// If the graph is disconnected, a minimum spanning forest is returned
// and `connected` is set to false.
inline MSTResult kruskal(const std::vector<std::vector<std::pair<int, int>>>& adj) {
    const int n = static_cast<int>(adj.size());

    // Gather all edges.
    std::vector<Edge> allEdges;
    for (int u = 0; u < n; ++u) {
        for (const auto& [v, w] : adj[u]) {
            allEdges.push_back({u, v, w});
        }
    }

    // Sort by weight ascending. Ties are broken by (smaller endpoint, larger endpoint)
    // so the edge order is a strict total order. With such an order the MST is unique,
    // which makes Kruskal and Prim (which uses the same order) pick identical edges.
    std::sort(allEdges.begin(), allEdges.end(), [](const Edge& a, const Edge& b) {
        if (a.w != b.w) return a.w < b.w;
        int aLo = std::min(a.u, a.v), bLo = std::min(b.u, b.v);
        if (aLo != bLo) return aLo < bLo;
        return std::max(a.u, a.v) < std::max(b.u, b.v);
    });

    MSTResult result;
    DisjointSet dsu(n);

    for (const Edge& e : allEdges) {
        if (dsu.unite(e.u, e.v)) {
            result.edges.push_back(e);
            result.totalWeight += e.w;
            if (static_cast<int>(result.edges.size()) == n - 1) break;
        }
    }

    result.connected = (n <= 1) || (dsu.componentCount() == 1);
    return result;
}
#endif // KRUSKAL_H
