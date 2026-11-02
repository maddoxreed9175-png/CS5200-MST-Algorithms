#ifndef PRIM_H
#define PRIM_H

#include <algorithm>
#include <functional>
#include <iostream>
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

// Shared types (also defined in kruskal.h; guarded so both headers can be included together).
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

// Prim's algorithm using an adjacency list and a min-priority queue (lazy version).
// adj[u] holds pairs (v, weight) for every edge u -- v; for an undirected graph
// each edge should appear in both adj[u] and adj[v].
// If the graph is disconnected, Prim is restarted from each unvisited vertex,
// producing a minimum spanning forest, and `connected` is set to false.
inline MSTResult prim(const std::vector<std::vector<std::pair<int, int>>>& adj,
                      int start = 0) {
    const int n = static_cast<int>(adj.size());
    MSTResult result;
    if (n == 0) return result;

    // Heap entry: (weight, to, from), ordered by smallest weight first.
    using Entry = std::tuple<int, int, int>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

    std::vector<bool> inTree(n, false);
    int components = 0;

    auto grow = [&](int src) {
        ++components;
        inTree[src] = true;
        for (const auto& [v, w] : adj[src]) {
            if (!inTree[v]) pq.emplace(w, v, src);
        }
        while (!pq.empty()) {
            auto [w, v, from] = pq.top();
            pq.pop();
            if (inTree[v]) continue;  // stale entry, vertex already reached cheaper
            inTree[v] = true;
            result.edges.push_back({from, v, w});
            result.totalWeight += w;
            for (const auto& [x, wx] : adj[v]) {
                if (!inTree[x]) pq.emplace(wx, x, v);
            }
        }
    };

    grow(start);
    for (int i = 0; i < n; ++i) {
        if (!inTree[i]) grow(i);
    }

    result.connected = (components == 1);
    return result;
}

#endif // PRIM_H
