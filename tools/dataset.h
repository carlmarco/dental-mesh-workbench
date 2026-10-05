#pragma once
// Readers for the local Teeth3DS+ / 3DTeethLand files (license-restricted, local only: D60).

#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "core/json.h"
#include "core/mesh.h"

namespace dmw::dataset {

inline std::string read_text(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// All files under `dirs` whose name ends with `suffix`, keyed by name with the suffix removed.
inline std::map<std::string, std::filesystem::path> index_files(const std::vector<std::string>& dirs,
                                                               const std::string& suffix) {
    std::map<std::string, std::filesystem::path> out;
    for (const auto& d : dirs) {
        for (const auto& e : std::filesystem::recursive_directory_iterator(d)) {
            const std::string name = e.path().filename().string();
            if (e.is_regular_file() && name.size() > suffix.size() &&
                name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
                out[name.substr(0, name.size() - suffix.size())] = e.path();
            }
        }
    }
    return out;
}

// 3DTeethLand landmarks of one class ("Cusp", "Mesial", ...) from a __kpt.json file.
inline std::vector<Vec3> read_landmarks(const std::filesystem::path& p, const std::string& cls) {
    std::vector<Vec3> out;
    const JsonResult r = parse_json(read_text(p));
    if (!r.ok()) return out;
    const JsonValue* objects = r.value.find("objects");
    if (!objects) return out;
    for (const JsonValue& o : objects->array) {
        const JsonValue* c = o.find("class");
        const JsonValue* xyz = o.find("coord");
        if (c && xyz && c->string == cls && xyz->array.size() == 3) {
            out.push_back({xyz->array[0].number, xyz->array[1].number, xyz->array[2].number});
        }
    }
    return out;
}

}  // namespace dmw::dataset
