#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace dmw {

struct JsonMember;

// Minimal JSON value (D65): enough for dataset annotations (Teeth3DS labels, 3DTeethLand landmarks).
// Text in, value out, no I/O (D14); errors are returned, not thrown (D15).
struct JsonValue {
    enum class Type { Null, Bool, Number, String, Array, Object };
    Type type = Type::Null;
    bool boolean = false;
    double number = 0.0;
    std::string string;
    std::vector<JsonValue> array;
    // Insertion order kept. A vector of a forward-declared struct (valid since C++17): a vector of
    // std::pair<std::string, JsonValue> instantiates the pair while JsonValue is still incomplete, which
    // libstdc++ rejects (libc++ happens to accept it; found by the first Linux CI build, D93).
    std::vector<JsonMember> object;

    // Member lookup; nullptr if absent or not an object.
    const JsonValue* find(std::string_view key) const;
};

struct JsonMember {
    std::string key;
    JsonValue value;
};

struct JsonResult {
    JsonValue value;
    std::string error;  // empty on success; otherwise "offset N: <reason>"
    bool ok() const { return error.empty(); }
};

JsonResult parse_json(std::string_view text);

}  // namespace dmw
