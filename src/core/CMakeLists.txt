add_library(dmw_core STATIC io_obj.cpp io_stl.cpp weld.cpp halfedge.cpp generate.cpp topology.cpp curvature.cpp sparse.cpp dec.cpp geodesic.cpp cholesky.cpp intrinsic.cpp io_write.cpp handles.cpp json.cpp metrics.cpp cusps.cpp)

# PUBLIC include dir: anything linking dmw_core can #include "core/<header>.h".
# The header lives under include/core/ so includes are namespaced by path.
target_include_directories(dmw_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Warnings are attached to the target (not set globally) so they don't leak
# into third-party code like Catch2, which we don't control.
target_compile_options(dmw_core PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
