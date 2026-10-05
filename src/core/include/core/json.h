#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace dmw {

// Minimal JSON value (D65): enough for dataset annotations (Teeth3DS labels, 3DTeethLand landmarks).
// Text in, value out, no I/O (D14); errors are returned, not thrown (D15).
struct JsonValue {
    enum class Type { Null, Bool, Number, String, Array, Object };
    Type type = Type::Null;
    bool boolean = false;
    double number = 0.0;
    std::string string;
    std::vector<JsonValue> array;
    std::vector<std::pair<std::string, JsonValue>> object;  // insertion order kept

    // Member lookup; nullptr if absent or not an object.
    const JsonValue* find(std::string_view key) const;
};

struct JsonResult {
    JsonValue value;
    std::string error;  // empty on success; otherwise "offset N: <reason>"
    bool ok() const { return error.empty(); }
};

JsonResult parse_json(std::string_view text);

}  // namespace dmw
