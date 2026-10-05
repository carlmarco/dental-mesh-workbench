add_executable(core_tests test_obj.cpp test_stl.cpp test_weld.cpp test_halfedge.cpp test_generate.cpp test_topology.cpp test_curvature.cpp test_sparse.cpp test_dec.cpp test_geodesic.cpp test_cholesky.cpp test_intrinsic.cpp test_io_write.cpp)

# Link against the core library plus Catch2 (without its own main()).
# Catch2WithMain supplies main(), so tests are just TEST_CASEs.
target_link_libraries(core_tests PRIVATE dmw_core Catch2::Catch2WithMain)

# catch_discover_tests asks the built test binary for its test cases and
# registers each one with CTest, so `ctest` reports them individually.
include(Catch)
catch_discover_tests(core_tests)
