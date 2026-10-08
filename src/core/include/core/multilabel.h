#pragma once

#include <cstdint>
#include <functional>
#include <span>
#include <vector>

namespace dmw {

// Multi-label energy minimization with a Potts pairwise term (D85):
//   E(f) = sum_i unary[i * labels + f_i] + sum_{edges (a, b, w)} w * [f_a != f_b],  w >= 0.
// alpha_expansion (Boykov, Veksler, Zabih, PAMI 23(11), 2001): for each label alpha in turn, solve the
// binary move "keep your label or switch to alpha" exactly by s-t min cut (Kolmogorov & Zabih 2004
// construction; the move energy is submodular for a metric such as Potts). E never increases; the result is
// a local minimum w.r.t. all expansion moves and, for Potts, within a factor 2 of the global optimum.
struct PottsEdge {
    std::uint32_t a, b;
    double w;
};

double potts_energy(std::uint32_t labels, std::span<const double> unary, std::span<const PottsEdge> edges,
                    std::span<const std::uint32_t> f);

// Optional candidate filter: fills mask[i] = 1 for nodes allowed to take `alpha` in this move (given the
// current labelling). Others keep their label. Nodes already labelled alpha are always left out of the move
// graph (exact: their choice does not matter); the filter is an approximation that shrinks the graph.
using ExpansionCandidates = std::function<void(std::uint32_t alpha, const std::vector<std::uint32_t>& f, std::vector<std::uint8_t>& mask)>;

// Improves `f` in place (initial labelling in, result out); stops after `max_sweeps` sweeps over all labels or
// when a sweep changes nothing. Returns the final energy.
double alpha_expansion(std::uint32_t labels, std::span<const double> unary, std::span<const PottsEdge> edges,
                       std::vector<std::uint32_t>& f, int max_sweeps, const ExpansionCandidates& candidates = {});

}  // namespace dmw
