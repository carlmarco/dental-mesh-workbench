# Benchmark tool (native only). Build it Release for meaningful numbers:
#   cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release -j
#   ./build-release/tools/dmw_bench build-release/bench
add_executable(dmw_bench bench.cpp)
target_link_libraries(dmw_bench PRIVATE dmw_core)
target_compile_options(dmw_bench PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

# Real-scan QA (reads local, license-restricted data; see D60, D62):
#   ./build-release/tools/scan_qa build-release/scan_qa.csv data/data_part_7
add_executable(scan_qa scan_qa.cpp)
target_link_libraries(scan_qa PRIVATE dmw_core)
target_compile_options(scan_qa PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

# Cusp detection evaluation (local, license-restricted data; D60, D67):
#   ./build-release/tools/cusp_eval sweep data/3DTeethLand_landmarks_train data/data_part_1
#   ./build-release/tools/cusp_eval test  data/3DTeethLand_landmarks_test  data/data_part_7
add_executable(cusp_eval cusp_eval.cpp)
target_link_libraries(cusp_eval PRIVATE dmw_core)
target_compile_options(cusp_eval PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

# Margin evaluation (local, license-restricted data; D60, D71-D73):
#   ./build-release/tools/margin_eval sweep data/data_part_1 5
#   ./build-release/tools/margin_eval test  data/data_part_5
add_executable(margin_eval margin_eval.cpp)
target_link_libraries(margin_eval PRIVATE dmw_core)
target_compile_options(margin_eval PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

# Margin failure analysis on VALIDATION scans (training scans not used in the sweeps; D74):
#   ./build-release/tools/margin_analyze data/data_part_1 build-release/margin_val.csv 5 0
add_executable(margin_analyze margin_analyze.cpp)
target_link_libraries(margin_analyze PRIVATE dmw_core)
target_compile_options(margin_analyze PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

# Cusp-seed classifier training (local data; model file stays local, D76):
#   ./build-release/tools/seed_train data/data_part_1 data/models/cusp_seed_logistic.json
add_executable(seed_train seed_train.cpp)
target_link_libraries(seed_train PRIVATE dmw_core)
target_compile_options(seed_train PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

# Learned per-vertex data term for the margin cut (D88): trains on Teeth3DS part 1 folds {1,2}.
#   ./build-release/tools/vertex_train data/data_part_1 data/models/vertex_logistic.json
add_executable(vertex_train vertex_train.cpp)
target_link_libraries(vertex_train PRIVATE dmw_core)
target_compile_options(vertex_train PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

# Undercut / path-of-insertion analysis on scans (M10b, D90).
#   ./build-release/tools/undercut_eval data/data_part_1 15
add_executable(undercut_eval undercut_eval.cpp)
target_link_libraries(undercut_eval PRIVATE dmw_core)
target_compile_options(undercut_eval PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
