#include "core/metrics.h"

#include <algorithm>
#include <cmath>
#include <tuple>

namespace dmw {

void MatchResult::add(const MatchResult& o) {
    true_positives += o.true_positives;
    detections += o.detections;
    ground_truth += o.ground_truth;
    matched_distances.insert(matched_distances.end(), o.matched_distances.begin(), o.matched_distances.end());
}

MatchResult match_points(std::span<const Vec3> det, std::span<const Vec3> gt, double tolerance) {
    MatchResult r;
    r.detections = det.size();
    r.ground_truth = gt.size();
    // All candidate pairs within tolerance, shortest first; each point used at most once.
    // O(|det| * |gt|): tens of points per scan, so brute force is fine.
    std::vector<std::tuple<double, std::size_t, std::size_t>> pairs;
    for (std::size_t i = 0; i < det.size(); ++i) {
        for (std::size_t j = 0; j < gt.size(); ++j) {
            const double dx = det[i].x - gt[j].x, dy = det[i].y - gt[j].y, dz = det[i].z - gt[j].z;
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (d <= tolerance) pairs.emplace_back(d, i, j);
        }
    }
    std::sort(pairs.begin(), pairs.end());
    std::vector<bool> det_used(det.size(), false), gt_used(gt.size(), false);
    for (const auto& [d, i, j] : pairs) {
        if (det_used[i] || gt_used[j]) continue;
        det_used[i] = gt_used[j] = true;
        ++r.true_positives;
        r.matched_distances.push_back(d);
    }
    return r;
}

}  // namespace dmw
