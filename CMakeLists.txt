cmake_minimum_required(VERSION 3.12)

project(Metagenomic-sample-analysis VERSION 1.0.0)

add_executable(MyProgram)

target_sources(MyProgram PRIVATE src/main.cpp)

add_subdirectory(external/bioparser)

target_link_libraries(MyProgram PRIVATE bioparser)