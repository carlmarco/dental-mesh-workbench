#pragma once

#include <cstdint>
#include <vector>

namespace dmw {

// s-t minimum cut / maximum flow on a sparse directed graph (D78), by the Boykov-Kolmogorov algorithm
// (Y. Boykov, V. Kolmogorov, "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy
// Minimization in Vision", IEEE PAMI 26(9), 2004). Two search trees, S from the source and T from the
// sink, grow until they touch (an augmenting path); after augmenting, saturated tree arcs orphan nodes,
// which are re-adopted (or freed) without restarting the search. Worst case is not polynomially better
// than Ford-Fulkerson, but on low-degree, grid-like graphs (meshes) it is fast in practice.
//
// Usage: add terminal capacities and edges, then solve() once. source_side(v) is the minimal source
// set: nodes reachable from the source in the final residual graph. Capacities must be >= 0 and finite
// (use a large finite value for hard constraints).
class MaxFlow {
public:
    explicit MaxFlow(std::uint32_t nodes);

    // Capacity source -> v and v -> sink. Calls accumulate.
    void add_terminal(std::uint32_t v, double source_cap, double sink_cap);
    // Capacity u -> v and v -> u.
    void add_edge(std::uint32_t u, std::uint32_t v, double cap, double reverse_cap);

    double solve();  // returns the maximum flow (= minimum cut value)
    bool source_side(std::uint32_t v) const { return nodes_[v].tree == kSource; }

private:
    static constexpr std::uint32_t kNone = 0xFFFFFFFFu;      // no arc / not in a tree
    static constexpr std::uint32_t kTerminal = 0xFFFFFFFEu;  // parent is the terminal itself
    static constexpr std::uint32_t kOrphan = 0xFFFFFFFDu;    // lost its parent arc, awaiting adoption
    static constexpr std::uint8_t kFree = 0, kSource = 1, kSink = 2;

    struct Arc {
        std::uint32_t head;  // node the arc points to
        std::uint32_t next;  // next arc leaving the same tail
        double residual;     // remaining capacity; the reverse arc is index ^ 1
    };
    struct Node {
        std::uint32_t first = kNone;   // first outgoing arc
        std::uint32_t parent = kNone;  // arc to the parent (node -> parent), kTerminal, kOrphan or kNone
        std::uint32_t timestamp = 0;   // when dist was last known exact (distance heuristic)
        std::uint32_t dist = 0;        // tree depth estimate
        double terminal = 0.0;         // residual: > 0 to the source, < 0 to the sink
        std::uint8_t tree = kFree;
        bool active = false;
    };

    void activate(std::uint32_t v);
    void augment(std::uint32_t middle);
    void adopt(std::uint32_t v);

    std::vector<Node> nodes_;
    std::vector<Arc> arcs_;
    std::vector<std::uint32_t> active_;  // FIFO, consumed from active_head_
    std::size_t active_head_ = 0;
    std::vector<std::uint32_t> orphans_;
    std::uint32_t time_ = 0;
    double flow_ = 0.0;
};

}  // namespace dmw
