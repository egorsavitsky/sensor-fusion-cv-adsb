cmake_minimum_required(VERSION 3.10)
project(sensor_fusion LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(fusion src/fusion.cpp src/geometry.cpp src/log_reader.cpp src/simulation.cpp)
target_include_directories(fusion PUBLIC include)
if(MSVC)
    target_compile_options(fusion PRIVATE /W4 /utf-8)
else()
    target_compile_options(fusion PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(fusion_core src/main.cpp)
target_link_libraries(fusion_core PRIVATE fusion)

include(CTest)
if(BUILD_TESTING)
    add_executable(fusion_tests tests/test_fusion.cpp)
    target_link_libraries(fusion_tests PRIVATE fusion)
    add_test(NAME fusion_tests COMMAND fusion_tests)
    add_test(NAME vertical_model COMMAND fusion_core
        --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/vertical.tsv
        --calibration-frames 0 --horizontal-sigma 0 --vertical-sigma 0
        --coverage 1 --dropout 0 --delay 0)
    set_tests_properties(vertical_model PROPERTIES PASS_REGULAR_EXPRESSION "Matches: 6, correct: 6, wrong: 0")
endif()
