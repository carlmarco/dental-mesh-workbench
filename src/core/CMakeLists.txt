add_library(dmw_core STATIC arith.cpp)

# PUBLIC include dir: anything linking dmw_core can #include "core/arith.h".
# The header lives under include/core/ so includes are namespaced by path.
target_include_directories(dmw_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Warnings are attached to the target (not set globally) so they don't leak
# into third-party code like Catch2, which we don't control.
target_compile_options(dmw_core PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
