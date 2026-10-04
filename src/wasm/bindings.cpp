// Thin embind layer: exposes dmw_core to JS. No logic lives here.
#include <emscripten/bind.h>

#include "core/arith.h"

// Expands to a static initializer that runs when the module loads and fills
// embind's registry (name -> function pointer + type info). The JS wrappers
// are generated from this registry at runtime.
EMSCRIPTEN_BINDINGS(dmw) {
    emscripten::function("add", &dmw::add);
}
