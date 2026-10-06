#include "core/maxflow.h"

#include <algorithm>
#include <limits>

namespace dmw {

MaxFlow::MaxFlow(std::uint32_t nodes) : nodes_(nodes) {}

void MaxFlow::add_terminal(std::uint32_t v, double source_cap, double sink_cap) {
    // Only the difference matters for the cut: min(source, sink) units go straight through v.
    const double through = std::min(source_cap, sink_cap);
    flow_ += through;
    nodes_[v].terminal += source_cap - sink_cap;
}

void MaxFlow::add_edge(std::uint32_t u, std::uint32_t v, double cap, double reverse_cap) {
    const auto a = static_cast<std::uint32_t>(arcs_.size());  // u -> v at a, v -> u at a ^ 1
    arcs_.push_back({v, nodes_[u].first, cap});
    arcs_.push_back({u, nodes_[v].first, reverse_cap});
    nodes_[u].first = a;
    nodes_[v].first = a + 1;
}

void MaxFlow::activate(std::uint32_t v) {
    if (nodes_[v].active) return;
    nodes_[v].active = true;
    active_.push_back(v);
}

double MaxFlow::solve() {
    for (std::uint32_t v = 0; v < nodes_.size(); ++v) {
        Node& n = nodes_[v];
        if (n.terminal == 0.0) continue;
        n.tree = n.terminal > 0.0 ? kSource : kSink;
        n.parent = kTerminal, n.timestamp = 0, n.dist = 1;
        activate(v);
    }
    std::uint32_t current = kNone;  // node whose arcs are being scanned (kept across an augmentation)
    for (;;) {
        // --- Growth: take the next active node still in a tree.
        if (current != kNone && nodes_[current].parent == kNone) current = kNone;  // freed by adoption
        if (current == kNone) {
            while (active_head_ < active_.size()) {
                const std::uint32_t v = active_[active_head_++];
                nodes_[v].active = false;
                if (nodes_[v].parent != kNone) {
                    current = v;
                    break;
                }
            }
            if (active_head_ == active_.size()) active_.clear(), active_head_ = 0;
            if (current == kNone) break;  // no active nodes: the trees cannot grow, flow is maximum
        }
        const std::uint32_t i = current;
        std::uint32_t middle = kNone;  // an arc from the S tree to the T tree
        for (std::uint32_t a = nodes_[i].first; a != kNone; a = arcs_[a].next) {
            // S grows along arcs i -> j with residual capacity; T along arcs j -> i (the reverse, a ^ 1).
            const bool s = nodes_[i].tree == kSource;
            if ((s ? arcs_[a].residual : arcs_[a ^ 1].residual) <= 0.0) continue;
            const std::uint32_t j = arcs_[a].head;
            Node& nj = nodes_[j];
            if (nj.parent == kNone) {
                nj.tree = nodes_[i].tree;
                nj.parent = a ^ 1;  // arc j -> i
                nj.timestamp = nodes_[i].timestamp, nj.dist = nodes_[i].dist + 1;
                activate(j);
            } else if (nj.tree != nodes_[i].tree) {
                middle = s ? a : (a ^ 1);
                break;
            } else if (nj.timestamp <= nodes_[i].timestamp && nj.dist > nodes_[i].dist) {
                // Distance heuristic: re-parent j onto the shorter path through i.
                nj.parent = a ^ 1, nj.timestamp = nodes_[i].timestamp, nj.dist = nodes_[i].dist + 1;
            }
        }
        ++time_;
        if (middle == kNone) {
            current = kNone;  // i is exhausted; it stays in its tree as a passive node
            continue;
        }
        augment(middle);  // `current` stays i: its remaining arcs are scanned next unless it is freed
        // --- Adoption: process orphans until none remain (new ones are appended as they appear).
        for (std::size_t k = 0; k < orphans_.size(); ++k) adopt(orphans_[k]);
        orphans_.clear();
    }
    return flow_;
}

void MaxFlow::augment(std::uint32_t middle) {
    // Bottleneck along source -> ... -> tail(middle) -> head(middle) -> ... -> sink.
    double b = arcs_[middle].residual;
    for (std::uint32_t v = arcs_[middle ^ 1].head;;) {  // S side: arcs parent -> v are parent[v] ^ 1
        const std::uint32_t a = nodes_[v].parent;
        if (a == kTerminal) {
            b = std::min(b, nodes_[v].terminal);
            break;
        }
        b = std::min(b, arcs_[a ^ 1].residual);
        v = arcs_[a].head;
    }
    for (std::uint32_t v = arcs_[middle].head;;) {  // T side: arcs v -> parent are parent[v]
        const std::uint32_t a = nodes_[v].parent;
        if (a == kTerminal) {
            b = std::min(b, -nodes_[v].terminal);
            break;
        }
        b = std::min(b, arcs_[a].residual);
        v = arcs_[a].head;
    }
    // Push b; every arc that saturates orphans the node below it.
    arcs_[middle].residual -= b;
    arcs_[middle ^ 1].residual += b;
    auto orphan = [this](std::uint32_t v) {
        nodes_[v].parent = kOrphan;
        orphans_.push_back(v);
    };
    for (std::uint32_t v = arcs_[middle ^ 1].head;;) {
        const std::uint32_t a = nodes_[v].parent;
        if (a == kTerminal) {
            nodes_[v].terminal -= b;
            if (nodes_[v].terminal <= 0.0) orphan(v);
            break;
        }
        arcs_[a ^ 1].residual -= b;
        arcs_[a].residual += b;
        const std::uint32_t up = arcs_[a].head;
        if (arcs_[a ^ 1].residual <= 0.0) orphan(v);
        v = up;
    }
    for (std::uint32_t v = arcs_[middle].head;;) {
        const std::uint32_t a = nodes_[v].parent;
        if (a == kTerminal) {
            nodes_[v].terminal += b;
            if (nodes_[v].terminal >= 0.0) orphan(v);
            break;
        }
        arcs_[a].residual -= b;
        arcs_[a ^ 1].residual += b;
        const std::uint32_t up = arcs_[a].head;
        if (arcs_[a].residual <= 0.0) orphan(v);
        v = up;
    }
    flow_ += b;
}

void MaxFlow::adopt(std::uint32_t i) {
    const std::uint8_t tree = nodes_[i].tree;
    // A neighbour j can be the new parent if the arc between them has capacity in the tree's direction
    // and j's own path still reaches the terminal (walk up; timestamps cache verified distances).
    std::uint32_t best_arc = kNone, best_dist = std::numeric_limits<std::uint32_t>::max();
    for (std::uint32_t a = nodes_[i].first; a != kNone; a = arcs_[a].next) {
        const std::uint32_t j = arcs_[a].head;
        if (nodes_[j].tree != tree || nodes_[j].parent == kNone) continue;
        if ((tree == kSource ? arcs_[a ^ 1].residual : arcs_[a].residual) <= 0.0) continue;
        std::uint32_t d = 0;
        bool reaches = false;
        for (std::uint32_t k = j;;) {
            if (nodes_[k].timestamp == time_) {
                d += nodes_[k].dist, reaches = true;
                break;
            }
            const std::uint32_t p = nodes_[k].parent;
            ++d;
            if (p == kTerminal) {
                nodes_[k].timestamp = time_, nodes_[k].dist = 1, reaches = true;
                break;
            }
            if (p == kOrphan) break;
            k = arcs_[p].head;
        }
        if (!reaches) continue;
        if (d < best_dist) best_arc = a, best_dist = d;
        for (std::uint32_t k = j; nodes_[k].timestamp != time_; k = arcs_[nodes_[k].parent].head) {
            nodes_[k].timestamp = time_, nodes_[k].dist = d--;  // cache exact distances on the path
        }
    }
    if (best_arc != kNone) {
        nodes_[i].parent = best_arc, nodes_[i].timestamp = time_, nodes_[i].dist = best_dist + 1;
        return;
    }
    // No parent: i becomes free. Neighbours that could reach it become active (they may regrow into
    // this region); its children become orphans.
    for (std::uint32_t a = nodes_[i].first; a != kNone; a = arcs_[a].next) {
        const std::uint32_t j = arcs_[a].head;
        if (nodes_[j].tree != tree || nodes_[j].parent == kNone) continue;
        if ((tree == kSource ? arcs_[a ^ 1].residual : arcs_[a].residual) > 0.0) activate(j);
        const std::uint32_t p = nodes_[j].parent;
        if (p != kTerminal && p != kOrphan && arcs_[p].head == i) {
            nodes_[j].parent = kOrphan;
            orphans_.push_back(j);
        }
    }
    nodes_[i].parent = kNone, nodes_[i].tree = kFree;
}

}  // namespace dmw
