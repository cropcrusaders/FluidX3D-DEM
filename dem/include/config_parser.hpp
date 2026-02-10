#pragma once
// JSON configuration parser for DEM simulator
// Uses a minimal JSON parser (included) to avoid nlohmann/json dependency

#include "simulator.hpp"
#include <string>

namespace dem {

// Minimal JSON value type (self-contained, no external deps)
class JsonValue {
public:
    enum Type { NUL, BOOL, NUMBER, STRING, ARRAY, OBJECT };

    JsonValue() : type_(NUL) {}
    JsonValue(bool b) : type_(BOOL), bool_val_(b) {}
    JsonValue(double d) : type_(NUMBER), num_val_(d) {}
    JsonValue(const std::string& s) : type_(STRING), str_val_(s) {}
    JsonValue(const char* s) : type_(STRING), str_val_(s) {}

    Type type() const { return type_; }
    bool is_null() const { return type_ == NUL; }
    bool is_bool() const { return type_ == BOOL; }
    bool is_number() const { return type_ == NUMBER; }
    bool is_string() const { return type_ == STRING; }
    bool is_array() const { return type_ == ARRAY; }
    bool is_object() const { return type_ == OBJECT; }

    bool as_bool(bool def = false) const { return type_ == BOOL ? bool_val_ : def; }
    double as_number(double def = 0.0) const { return type_ == NUMBER ? num_val_ : def; }
    int as_int(int def = 0) const { return type_ == NUMBER ? (int)num_val_ : def; }
    const std::string& as_string(const std::string& def = "") const {
        static std::string empty;
        return type_ == STRING ? str_val_ : (def.empty() ? empty : def);
    }

    // Array access
    size_t size() const { return arr_val_.size(); }
    const JsonValue& operator[](size_t i) const {
        static JsonValue null_val;
        return i < arr_val_.size() ? arr_val_[i] : null_val;
    }

    // Object access
    const JsonValue& operator[](const std::string& key) const {
        static JsonValue null_val;
        for (auto& p : obj_val_)
            if (p.first == key) return p.second;
        return null_val;
    }
    bool has(const std::string& key) const {
        for (auto& p : obj_val_)
            if (p.first == key) return true;
        return false;
    }

    // Builder methods
    void push_back(const JsonValue& v) { type_ = ARRAY; arr_val_.push_back(v); }
    void set(const std::string& key, const JsonValue& v) {
        type_ = OBJECT;
        for (auto& p : obj_val_) {
            if (p.first == key) { p.second = v; return; }
        }
        obj_val_.push_back({key, v});
    }

    const std::vector<std::pair<std::string, JsonValue>>& object_items() const { return obj_val_; }
    const std::vector<JsonValue>& array_items() const { return arr_val_; }

    // Parse from string
    static JsonValue parse(const std::string& json);

private:
    Type type_ = NUL;
    bool bool_val_ = false;
    double num_val_ = 0.0;
    std::string str_val_;
    std::vector<JsonValue> arr_val_;
    std::vector<std::pair<std::string, JsonValue>> obj_val_;
};

// Parse a JSON config file and populate simulator
bool load_config(const std::string& path, Simulator& sim);

// Parse JSON string
bool parse_config(const std::string& json_str, Simulator& sim);

} // namespace dem
