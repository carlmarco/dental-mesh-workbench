#include "core/json.h"

#include <cstdint>
#include <string>

#include "detail/parse_number.h"

namespace dmw {

const JsonValue* JsonValue::find(std::string_view key) const {
    if (type != Type::Object) return nullptr;
    for (const auto& [k, v] : object) {
        if (k == key) return &v;
    }
    return nullptr;
}

namespace {

// Recursive descent over the JSON grammar. Depth-limited so hostile input can't overflow the stack.
class Parser {
public:
    explicit Parser(std::string_view text) : s_(text) {}

    JsonResult run() {
        JsonResult r;
        skip_ws();
        if (!value(r.value, 0)) {
            r.error = "offset " + std::to_string(pos_) + ": " + error_;
            r.value = {};
            return r;
        }
        skip_ws();
        if (pos_ != s_.size()) r.error = "offset " + std::to_string(pos_) + ": trailing characters";
        return r;
    }

private:
    static constexpr int kMaxDepth = 256;

    bool fail(const char* why) {
        error_ = why;
        return false;
    }
    void skip_ws() {
        while (pos_ < s_.size() && (s_[pos_] == ' ' || s_[pos_] == '\t' || s_[pos_] == '\n' || s_[pos_] == '\r')) ++pos_;
    }
    bool literal(std::string_view word) {
        if (s_.substr(pos_, word.size()) != word) return fail("invalid literal");
        pos_ += word.size();
        return true;
    }

    bool value(JsonValue& out, int depth) {
        if (depth > kMaxDepth) return fail("nesting too deep");
        if (pos_ >= s_.size()) return fail("unexpected end of input");
        const char c = s_[pos_];
        if (c == '{') return object(out, depth);
        if (c == '[') return array(out, depth);
        if (c == '"') {
            out.type = JsonValue::Type::String;
            return string(out.string);
        }
        if (c == 't') return out.type = JsonValue::Type::Bool, out.boolean = true, literal("true");
        if (c == 'f') return out.type = JsonValue::Type::Bool, out.boolean = false, literal("false");
        if (c == 'n') return out.type = JsonValue::Type::Null, literal("null");
        return number(out);
    }

    bool number(JsonValue& out) {
        const std::size_t begin = pos_;
        while (pos_ < s_.size() && (std::string_view("+-0123456789.eE").find(s_[pos_]) != std::string_view::npos)) ++pos_;
        if (!detail::parse_double(s_.substr(begin, pos_ - begin), out.number)) return fail("invalid number");
        out.type = JsonValue::Type::Number;
        return true;
    }

    bool string(std::string& out) {
        ++pos_;  // opening quote
        while (pos_ < s_.size() && s_[pos_] != '"') {
            char c = s_[pos_++];
            if (c == '\\') {
                if (pos_ >= s_.size()) return fail("unterminated escape");
                const char e = s_[pos_++];
                switch (e) {
                    case '"': c = '"'; break;
                    case '\\': c = '\\'; break;
                    case '/': c = '/'; break;
                    case 'b': c = '\b'; break;
                    case 'f': c = '\f'; break;
                    case 'n': c = '\n'; break;
                    case 'r': c = '\r'; break;
                    case 't': c = '\t'; break;
                    case 'u': {  // \uXXXX: keep ASCII, replace anything else (not needed for our data)
                        if (pos_ + 4 > s_.size()) return fail("short \\u escape");
                        std::uint32_t code = 0;
                        for (int i = 0; i < 4; ++i) {
                            const char h = s_[pos_++];
                            code <<= 4;
                            if (h >= '0' && h <= '9') code |= static_cast<std::uint32_t>(h - '0');
                            else if (h >= 'a' && h <= 'f') code |= static_cast<std::uint32_t>(h - 'a' + 10);
                            else if (h >= 'A' && h <= 'F') code |= static_cast<std::uint32_t>(h - 'A' + 10);
                            else return fail("bad \\u escape");
                        }
                        c = code < 0x80 ? static_cast<char>(code) : '?';
                        break;
                    }
                    default: return fail("unknown escape");
                }
            }
            out.push_back(c);
        }
        if (pos_ >= s_.size()) return fail("unterminated string");
        ++pos_;  // closing quote
        return true;
    }

    bool array(JsonValue& out, int depth) {
        out.type = JsonValue::Type::Array;
        ++pos_;
        skip_ws();
        if (pos_ < s_.size() && s_[pos_] == ']') return ++pos_, true;
        for (;;) {
            out.array.emplace_back();
            skip_ws();
            if (!value(out.array.back(), depth + 1)) return false;
            skip_ws();
            if (pos_ < s_.size() && s_[pos_] == ',') { ++pos_; continue; }
            if (pos_ < s_.size() && s_[pos_] == ']') return ++pos_, true;
            return fail("expected ',' or ']'");
        }
    }

    bool object(JsonValue& out, int depth) {
        out.type = JsonValue::Type::Object;
        ++pos_;
        skip_ws();
        if (pos_ < s_.size() && s_[pos_] == '}') return ++pos_, true;
        for (;;) {
            skip_ws();
            if (pos_ >= s_.size() || s_[pos_] != '"') return fail("expected a string key");
            std::string key;
            if (!string(key)) return false;
            skip_ws();
            if (pos_ >= s_.size() || s_[pos_] != ':') return fail("expected ':'");
            ++pos_;
            skip_ws();
            out.object.emplace_back(std::move(key), JsonValue{});
            if (!value(out.object.back().second, depth + 1)) return false;
            skip_ws();
            if (pos_ < s_.size() && s_[pos_] == ',') { ++pos_; continue; }
            if (pos_ < s_.size() && s_[pos_] == '}') return ++pos_, true;
            return fail("expected ',' or '}'");
        }
    }

    std::string_view s_;
    std::size_t pos_ = 0;
    std::string error_;
};

}  // namespace

JsonResult parse_json(std::string_view text) { return Parser(text).run(); }

}  // namespace dmw
