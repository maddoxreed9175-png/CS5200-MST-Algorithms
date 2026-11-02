#ifndef TEST_H
#define TEST_H

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "generate.h"  // provides generateConnectedGraph(V, E)
#include "kruskal.h"
#include "prim.h"

// Testing procedure for kruskal() and prim().
//
// For every test graph the harness:
//   1. runs Kruskal and Prim,
//   2. checks each result independently: exactly V-1 edges, every edge really
//      exists in the graph with that weight, no cycles, all vertices connected,
//      and totalWeight equals the sum of the edge weights,
//   3. compares the two total weights (and an expected weight where one is known,
//      and a brute-force optimum for small graphs),
//   4. writes printMST() output of Kruskal to one file and of Prim to another.
// At the end the two output files are compared line by line.
//
// Usage:  bool ok = mst_test::runAllTests();

namespace mst_test {

using Graph = std::vector<std::vector<std::pair<int, int>>>;

struct Report {
    int cases = 0;
    int failedCases = 0;
    std::vector<std::string> mismatches;
};

inline void addEdge(Graph& g, int u, int v, int w) {
    g[u].push_back({v, w});
    g[v].push_back({u, w});
}

// Independent verification that `r` is a spanning tree of `g`.
inline bool validateTree(const Graph& g, const MSTResult& r, std::string& err) {
    const int n = static_cast<int>(g.size());
    const size_t expectedEdges = (n == 0) ? 0 : static_cast<size_t>(n - 1);

    if (r.edges.size() != expectedEdges) {
        err = "expected " + std::to_string(expectedEdges) + " edges, got " +
              std::to_string(r.edges.size());
        return false;
    }

    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    auto find = [&](int x) {
        while (p[x] != x) {
            p[x] = p[p[x]];
            x = p[x];
        }
        return x;
    };

    long long sum = 0;
    for (const Edge& e : r.edges) {
        if (e.u < 0 || e.u >= n || e.v < 0 || e.v >= n) {
            err = "edge has vertex out of range";
            return false;
        }
        bool exists = false;
        for (const auto& [v, w] : g[e.u]) exists = exists || (v == e.v && w == e.w);
        for (const auto& [v, w] : g[e.v]) exists = exists || (v == e.u && w == e.w);
        if (!exists) {
            err = "edge " + std::to_string(e.u) + "-" + std::to_string(e.v) + " (w=" +
                  std::to_string(e.w) + ") is not in the graph";
            return false;
        }
        int ru = find(e.u), rv = find(e.v);
        if (ru == rv) {
            err = "edge " + std::to_string(e.u) + "-" + std::to_string(e.v) + " creates a cycle";
            return false;
        }
        p[ru] = rv;
        sum += e.w;
    }

    for (int i = 1; i < n; ++i) {
        if (find(i) != find(0)) {
            err = "vertex " + std::to_string(i) + " is not connected (not spanning)";
            return false;
        }
    }
    if (sum != r.totalWeight) {
        err = "totalWeight " + std::to_string(r.totalWeight) + " != sum of edges " +
              std::to_string(sum);
        return false;
    }
    if (!r.connected) {
        err = "result flagged as disconnected for a connected graph";
        return false;
    }
    return true;
}

// Exhaustive MST weight for small graphs (tries every subset of V-1 edges).
inline long long bruteForceMST(const Graph& g) {
    const int n = static_cast<int>(g.size());
    if (n <= 1) return 0;

    std::set<std::tuple<int, int, int>> uniq;
    for (int u = 0; u < n; ++u)
        for (const auto& [v, w] : g[u])
            if (u != v) uniq.insert({std::min(u, v), std::max(u, v), w});
    std::vector<std::tuple<int, int, int>> es(uniq.begin(), uniq.end());
    const int m = static_cast<int>(es.size());

    long long best = -1;
    for (uint32_t mask = 0; mask < (1u << m); ++mask) {
        if (static_cast<int>(std::bitset<32>(mask).count()) != n - 1) continue;
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        auto find = [&](int x) {
            while (p[x] != x) x = p[x];
            return x;
        };
        long long sum = 0;
        bool ok = true;
        for (int i = 0; i < m && ok; ++i) {
            if (!((mask >> i) & 1u)) continue;
            auto [u, v, w] = es[i];
            int ru = find(u), rv = find(v);
            if (ru == rv) ok = false;
            else { p[ru] = rv; sum += w; }
        }
        if (ok && (best < 0 || sum < best)) best = sum;
    }
    return best;
}

// Runs both algorithms on one graph, validates, compares, and logs printed output.
// expected < 0 means "no hand-computed weight available".
inline void runCase(const std::string& name, const Graph& g, long long expected,
                    bool useBruteForce, std::ostream& kOut, std::ostream& pOut,
                    Report& rep) {
    ++rep.cases;
    MSTResult k = kruskal(g);
    MSTResult p = prim(g);

    kOut << "=== " << name << " ===\n";
    printMST(k, kOut);
    pOut << "=== " << name << " ===\n";
    printMST(p, pOut);

    std::vector<std::string> problems;
    std::string err;
    if (!validateTree(g, k, err)) problems.push_back("Kruskal result invalid: " + err);
    if (!validateTree(g, p, err)) problems.push_back("Prim result invalid: " + err);

    if (k.totalWeight != p.totalWeight) {
        problems.push_back("weight mismatch: Kruskal=" + std::to_string(k.totalWeight) +
                           " Prim=" + std::to_string(p.totalWeight));
    }
    if (expected >= 0) {
        if (k.totalWeight != expected)
            problems.push_back("Kruskal weight " + std::to_string(k.totalWeight) +
                               " != expected " + std::to_string(expected));
        if (p.totalWeight != expected)
            problems.push_back("Prim weight " + std::to_string(p.totalWeight) +
                               " != expected " + std::to_string(expected));
    }
    if (useBruteForce) {
        long long b = bruteForceMST(g);
        if (k.totalWeight != b)
            problems.push_back("Kruskal weight " + std::to_string(k.totalWeight) +
                               " != brute-force optimum " + std::to_string(b));
        if (p.totalWeight != b)
            problems.push_back("Prim weight " + std::to_string(p.totalWeight) +
                               " != brute-force optimum " + std::to_string(b));
    }

    if (!problems.empty()) {
        ++rep.failedCases;
        for (const auto& pr : problems) rep.mismatches.push_back(name + ": " + pr);
    }
}

// Line-by-line comparison of two files. Returns the number of differing lines.
inline int compareFiles(const std::string& fileA, const std::string& fileB) {
    std::ifstream a(fileA), b(fileB);
    std::string la, lb;
    int line = 0, diffs = 0;
    while (true) {
        bool ha = static_cast<bool>(std::getline(a, la));
        bool hb = static_cast<bool>(std::getline(b, lb));
        if (!ha && !hb) break;
        ++line;
        if (!ha) la = "<end of file>";
        if (!hb) lb = "<end of file>";
        if (la != lb) {
            if (diffs < 5) {
                std::cout << "  line " << line << ": " << fileA << " has \"" << la << "\", "
                          << fileB << " has \"" << lb << "\"\n";
            }
            ++diffs;
        }
    }
    return diffs;
}

inline bool runAllTests(const std::string& kruskalFile = "kruskal_output.txt",
                        const std::string& primFile = "prim_output.txt") {
    std::ofstream kOut(kruskalFile), pOut(primFile);
    Report rep;

    // ---------- 1. Small graphs with manually verifiable MSTs ----------
    {
        Graph g(1);
        runCase("manual: single vertex (MST weight 0)", g, 0, true, kOut, pOut, rep);
    }
    {
        Graph g(2);
        addEdge(g, 0, 1, 7);
        runCase("manual: two vertices (MST weight 7)", g, 7, true, kOut, pOut, rep);
    }
    {   // Triangle: pick 1 and 2, skip 3.
        Graph g(3);
        addEdge(g, 0, 1, 1); addEdge(g, 1, 2, 2); addEdge(g, 0, 2, 3);
        runCase("manual: triangle (MST weight 3)", g, 3, true, kOut, pOut, rep);
    }
    {   // Square with a heavy diagonal: 1+2+3 = 6.
        Graph g(4);
        addEdge(g, 0, 1, 1); addEdge(g, 1, 2, 2); addEdge(g, 2, 3, 3);
        addEdge(g, 3, 0, 4); addEdge(g, 0, 2, 5);
        runCase("manual: square + diagonal (MST weight 6)", g, 6, true, kOut, pOut, rep);
    }
    {   // Classic 5-vertex example: 2+3+5+6 = 16.
        Graph g(5);
        addEdge(g, 0, 1, 2); addEdge(g, 0, 3, 6); addEdge(g, 1, 2, 3);
        addEdge(g, 1, 3, 8); addEdge(g, 1, 4, 5); addEdge(g, 2, 4, 7);
        addEdge(g, 3, 4, 9);
        runCase("manual: 5-vertex classic (MST weight 16)", g, 16, true, kOut, pOut, rep);
    }
    {   // Already a tree: MST is the whole graph, 1+2+3+4 = 10.
        Graph g(5);
        addEdge(g, 0, 1, 1); addEdge(g, 1, 2, 2); addEdge(g, 2, 3, 3); addEdge(g, 3, 4, 4);
        runCase("manual: path graph (MST weight 10)", g, 10, true, kOut, pOut, rep);
    }

    // ---------- 2. Graphs containing equal-weight edges ----------
    {   // K4, every weight 1: any spanning tree has weight 3.
        Graph g(4);
        for (int u = 0; u < 4; ++u)
            for (int v = u + 1; v < 4; ++v) addEdge(g, u, v, 1);
        runCase("ties: K4 all weight 1 (MST weight 3)", g, 3, true, kOut, pOut, rep);
    }
    {   // K5, every weight 7: 4 * 7 = 28.
        Graph g(5);
        for (int u = 0; u < 5; ++u)
            for (int v = u + 1; v < 5; ++v) addEdge(g, u, v, 7);
        runCase("ties: K5 all weight 7 (MST weight 28)", g, 28, true, kOut, pOut, rep);
    }
    {   // Triangle of weight-1 edges plus two weight-2 edges: 1+1+2 = 4.
        Graph g(4);
        addEdge(g, 0, 1, 1); addEdge(g, 1, 2, 1); addEdge(g, 0, 2, 1);
        addEdge(g, 2, 3, 2); addEdge(g, 1, 3, 2);
        runCase("ties: mixed equal weights (MST weight 4)", g, 4, true, kOut, pOut, rep);
    }
    {   // 4-cycle plus chord, all weight 2: 3 * 2 = 6.
        Graph g(4);
        addEdge(g, 0, 1, 2); addEdge(g, 1, 2, 2); addEdge(g, 2, 3, 2);
        addEdge(g, 3, 0, 2); addEdge(g, 0, 2, 2);
        runCase("ties: 4-cycle + chord all weight 2 (MST weight 6)", g, 6, true, kOut, pOut, rep);
    }

    // ---------- 3. Sparse and dense connected graphs ----------
    {
        struct Spec { const char* label; int V, E; };
        const Spec specs[] = {
            {"sparse", 50, 55},   {"sparse", 200, 210}, {"sparse", 500, 520},
            {"dense", 30, 400},   {"dense", 50, 1000},  {"dense", 60, 1700},
        };
        for (const auto& s : specs) {
            Graph g = generateConnectedGraph(s.V, s.E);
            runCase(std::string(s.label) + ": V=" + std::to_string(s.V) +
                        " E=" + std::to_string(s.E),
                    g, -1, false, kOut, pOut, rep);
        }
    }

    // ---------- 4. Small random graphs checked against brute force ----------
    std::mt19937 rng(12345);
    for (int i = 1; i <= 30; ++i) {
        int V = 2 + static_cast<int>(rng() % 6);            // 2..7
        int maxE = std::min(V * (V - 1) / 2, 14);
        int E = (V - 1) + static_cast<int>(rng() % (maxE - (V - 1) + 1));
        Graph g = generateConnectedGraph(V, E);
        runCase("small-random #" + std::to_string(i) + " V=" + std::to_string(V) +
                    " E=" + std::to_string(E),
                g, -1, true, kOut, pOut, rep);
    }

    // ---------- 5. At least 100 random connected graphs ----------
    for (int i = 1; i <= 120; ++i) {
        int V = 2 + static_cast<int>(rng() % 79);           // 2..80
        int maxE = V * (V - 1) / 2;
        int E = (V - 1) + static_cast<int>(rng() % (maxE - (V - 1) + 1));
        Graph g = generateConnectedGraph(V, E);
        runCase("random #" + std::to_string(i) + " V=" + std::to_string(V) +
                    " E=" + std::to_string(E),
                g, -1, false, kOut, pOut, rep);
    }

    kOut.close();
    pOut.close();

    // ---------- Report ----------
    std::cout << "Test cases run: " << rep.cases << "\n";
    std::cout << "Cases with problems: " << rep.failedCases << "\n";
    for (const auto& m : rep.mismatches) std::cout << "  MISMATCH: " << m << "\n";

    std::cout << "Comparing " << kruskalFile << " and " << primFile << "...\n";
    int diffs = compareFiles(kruskalFile, primFile);
    if (diffs == 0) std::cout << "  Output files are identical.\n";
    else std::cout << "  Output files differ on " << diffs << " line(s).\n";

    bool ok = (rep.failedCases == 0 && diffs == 0);
    std::cout << (ok ? "ALL TESTS PASSED\n" : "SOME TESTS FAILED\n");
    return ok;
}

}  // namespace mst_test

#endif  // TEST_H
