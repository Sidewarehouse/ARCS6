cmake_minimum_required(VERSION 3.16)
project(ARCS_LIB)

include(${CMAKE_CURRENT_SOURCE_DIR}/../sys/ARCS_OPTIONS.cmake)

# ライブラリのソースを収集
file(GLOB_RECURSE ARCS_LIB CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.hh"
)

add_library(ARCS_LIB STATIC ${ARCS_LIB})
target_include_directories(ARCS_LIB
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../sys # sysに依存
)
# C++17機能を要求
target_compile_features(ARCS_LIB PUBLIC cxx_std_17)
# コンパイラオプションを設定
target_compile_options(ARCS_LIB PRIVATE ${ARCS_COMMON_CXX_OPTIONS})