#ifndef EVAL_H
#define EVAL_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "generate.h"  // provides generateConnectedGraph(V, E)
#include "kruskal.h"
#include "prim.h"

// Benchmark of kruskal() and prim() on graphs of increasing size and two densities.
//
// Usage:
//     #include "eval.h"
//     int main() { mst_eval::runBenchmarks(); }          // table to std::cout
//     // or: mst_eval::runBenchmarks(std::cout, "bench.csv");   // also writes a CSV
//
// Compile with optimizations (e.g. -O2), otherwise timings are not meaningful.
//
// Only the algorithm call is timed; graph generation is excluded. Each timing is the
// median of several repetitions.
//
// Theoretical complexity of the implementations in this project (V vertices, E edges):
//   Kruskal : O(E log E) = O(E log V)
//             sorting the edges dominates; the disjoint set with path compression and
//             union by size adds O(E * alpha(V)), where alpha is the inverse Ackermann
//             function (effectively constant).
//   Prim    : O(E log E) = O(E log V)
//             lazy binary-heap version: up to E pushes and pops at O(log E) each.
//             (An O(V^2) array-based Prim would be asymptotically better on very dense
//             graphs, but the task asks for a priority-queue implementation.)

namespace mst_eval {

using Graph = std::vector<std::vector<std::pair<int, int>>>;

struct Config {
    std::string density;  // "sparse" or "dense"
    int V;
    long long E;          // requested number of (undirected) edges
};

struct Row {
    std::string density;
    int V = 0;
    long long E = 0;        // actual undirected edge count (adjacency entries / 2)
    double kruskalMs = 0;
    double primMs = 0;
    long long kruskalWeight = 0;
    long long primWeight = 0;
};

template <class F>
double medianTimeMs(F algorithm, int reps, long long& sink) {
    std::vector<double> times;
    for (int r = 0; r < reps; ++r) {
        auto start = std::chrono::steady_clock::now();
        MSTResult res = algorithm();
        auto stop = std::chrono::steady_clock::now();
        sink += res.totalWeight + static_cast<long long>(res.edges.size());  // keep result alive
        times.push_back(std::chrono::duration<double, std::milli>(stop - start).count());
    }
    std::sort(times.begin(), times.end());
    return times[times.size() / 2];
}

// Default experiment: 4 sizes x 2 densities.
//   sparse: E = 4V          (average degree 8)
//   dense : E = V(V-1)/4    (half of all possible edges)
inline std::vector<Config> defaultConfigs() {
    std::vector<Config> cfgs;
    for (int V : {5000, 10000, 20000, 40000}) cfgs.push_back({"sparse", V, 4LL * V});
    for (int V : {250, 500, 1000, 2000})
        cfgs.push_back({"dense", V, static_cast<long long>(V) * (V - 1) / 4});
    return cfgs;
}

inline std::vector<Row> runExperiments(const std::vector<Config>& cfgs, int reps = 5) {
    std::vector<Row> rows;
    long long sink = 0;
    for (const Config& c : cfgs) {
        Graph g = generateConnectedGraph(c.V, static_cast<int>(c.E));

        long long entries = 0;
        for (const auto& nbrs : g) entries += static_cast<long long>(nbrs.size());

        Row row;
        row.density = c.density;
        row.V = c.V;
        row.E = entries / 2;
        row.kruskalMs = medianTimeMs([&] { return kruskal(g); }, reps, sink);
        row.primMs = medianTimeMs([&] { return prim(g); }, reps, sink);
        row.kruskalWeight = kruskal(g).totalWeight;
        row.primWeight = prim(g).totalWeight;
        rows.push_back(row);
    }
    if (sink == -1) std::cerr << "";  // never true; prevents the sink being optimized away
    return rows;
}

inline void printComplexity(std::ostream& out) {
    out << "Theoretical time complexity\n"
        << "  Kruskal (sort + DSU w/ path compression & union by size): "
           "O(E log E) = O(E log V)\n"
        << "  Prim    (lazy binary min-heap, adjacency list)           : "
           "O(E log E) = O(E log V)\n"
        << "  Both are O(E log V); the DSU's alpha(V) term is effectively constant.\n\n";
}

// "Model" column: E*log2(E) relative to the first row of the same density, i.e. how
// much slower the run should be than the smallest run if the O(E log E) model is exact.
inline void printReport(const std::vector<Row>& rows, std::ostream& out = std::cout) {
    printComplexity(out);

    out << "Runtimes are medians in milliseconds. 'x' columns are growth relative to the\n"
           "smallest graph of the same density; 'Model x' is the predicted growth E*log2(E).\n\n";

    out << std::left << std::setw(8) << "Density" << std::right << std::setw(8) << "V"
        << std::setw(11) << "E" << std::setw(13) << "Kruskal ms" << std::setw(11) << "Prim ms"
        << std::setw(10) << "K/P" << std::setw(11) << "Kruskal x" << std::setw(9) << "Prim x"
        << std::setw(10) << "Model x" << std::setw(8) << "Match" << "\n";
    out << std::string(99, '-') << "\n";

    const Row* base = nullptr;
    for (const Row& r : rows) {
        if (!base || base->density != r.density) base = &r;
        double model = (static_cast<double>(r.E) * std::log2(static_cast<double>(r.E))) /
                       (static_cast<double>(base->E) * std::log2(static_cast<double>(base->E)));
        bool match = (r.kruskalWeight == r.primWeight);

        out << std::left << std::setw(8) << r.density << std::right << std::setw(8) << r.V
            << std::setw(11) << r.E << std::fixed << std::setprecision(3) << std::setw(13)
            << r.kruskalMs << std::setw(11) << r.primMs << std::setprecision(2) << std::setw(10)
            << (r.primMs > 0 ? r.kruskalMs / r.primMs : 0.0) << std::setw(11)
            << (base->kruskalMs > 0 ? r.kruskalMs / base->kruskalMs : 0.0) << std::setw(9)
            << (base->primMs > 0 ? r.primMs / base->primMs : 0.0) << std::setw(10) << model
            << std::setw(8) << (match ? "yes" : "NO") << "\n";
    }
    out << "\nK/P > 1 means Prim was faster; K/P < 1 means Kruskal was faster.\n";
    out << "'Match' confirms both algorithms returned the same total MST weight.\n";
}

inline void writeCSV(const std::vector<Row>& rows, const std::string& path) {
    std::ofstream f(path);
    f << "density,V,E,kruskal_ms,prim_ms,kruskal_weight,prim_weight\n";
    for (const Row& r : rows) {
        f << r.density << ',' << r.V << ',' << r.E << ',' << r.kruskalMs << ',' << r.primMs << ','
          << r.kruskalWeight << ',' << r.primWeight << '\n';
    }
}

// Runs the default experiment, prints the report, optionally writes a CSV.
inline std::vector<Row> runBenchmarks(std::ostream& out = std::cout,
                                      const std::string& csvPath = "") {
    std::vector<Row> rows = runExperiments(defaultConfigs());
    printReport(rows, out);
    if (!csvPath.empty()) writeCSV(rows, csvPath);
    return rows;
}

}  // namespace mst_eval

#endif  // EVAL_H
