cmake_minimum_required(VERSION 3.12)

project(Metagenomic-sample-analysis VERSION 1.1.0)

add_subdirectory(external/bioparser)

set(SRC_FILES
    src/main.cpp
    src/loading_sequences.cpp
    src/minimizers.cpp
)

add_executable(distribution_analysis ${SRC_FILES})

target_include_directories(distribution_analysis PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    include
)

target_link_libraries(distribution_analysis PRIVATE bioparser)