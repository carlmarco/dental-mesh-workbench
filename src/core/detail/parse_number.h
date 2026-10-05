#pragma once
// Private to dmw_core (not under include/): shared by the text parsers.

#include <cmath>
#include <cstdlib>
#include <string>
#include <string_view>

namespace dmw::detail {

// std::from_chars<double> is unavailable in Apple's libc++ (DECISIONS.md D18), so use strtod.
// strtod needs a NUL-terminated string, and a string_view token points into the middle of
// the file, so it is copied first (short tokens fit in std::string's small buffer, no heap).
// Rejects partial parses ("1abc") and non-finite values (inf, nan, overflow).
inline bool parse_double(std::string_view tok, double& out) {
    if (tok.empty()) return false;
    const std::string buf(tok);
    char* end = nullptr;
    const double v = std::strtod(buf.c_str(), &end);
    if (end != buf.c_str() + buf.size() || !std::isfinite(v)) return false;
    out = v;
    return true;
}

}  // namespace dmw::detail
