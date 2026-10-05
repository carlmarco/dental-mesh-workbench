# Benchmark tool (native only). Build it Release for meaningful numbers:
#   cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release -j
#   ./build-release/tools/dmw_bench build-release/bench
add_executable(dmw_bench bench.cpp)
target_link_libraries(dmw_bench PRIVATE dmw_core)
target_compile_options(dmw_bench PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
