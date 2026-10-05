#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "core/mesh.h"

namespace dmw {

// Point-detection metrics (D66): one-to-one greedy matching of detections to ground truth by
// ascending distance, accepting pairs closer than `tolerance`.
struct MatchResult {
    std::size_t true_positives = 0;
    std::size_t detections = 0;
    std::size_t ground_truth = 0;
    std::vector<double> matched_distances;  // one per true positive

    double precision() const { return detections ? double(true_positives) / double(detections) : 0.0; }
    double recall() const { return ground_truth ? double(true_positives) / double(ground_truth) : 0.0; }
    double f1() const {
        const double p = precision(), r = recall();
        return p + r > 0.0 ? 2.0 * p * r / (p + r) : 0.0;
    }
    // Pool counts and distances across scans (micro-average).
    void add(const MatchResult& other);
};

MatchResult match_points(std::span<const Vec3> detections, std::span<const Vec3> ground_truth, double tolerance);

}  // namespace dmw
