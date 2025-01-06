cmake_minimum_required(VERSION 3.16)
project(ARCS)

set(CMAKE_CXX_STANDARD 17)
# 共通のコンパイルオプション
set(CMAKE_CXX_FLAGS "-fconcepts -Wall -Weffc++ -ftree-vectorize -march=native -DARCS_IN")
# Releaseのコンパイルオプション
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
# Debugのコンパイルオプション
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -Og")
# アセンブリリスト出力用
set(CMAKE_CXX_FLAGS_ASM  "-S -g")

# pthreadを見つける
find_package(Threads REQUIRED)

# ARCSのルートディレクトリを取得
get_filename_component(arcs_root_dir ${ARCS_SOURCE_DIR}/../../.. ABSOLUTE)
# lib, src, sysのファイルを.cmakeから取得
include(${arcs_root_dir}/lib/ARCS_LIB.cmake)
include(${arcs_root_dir}/sys/ARCS_SYS.cmake)

# ロボットフォルダを取得
get_filename_component(arcs_robot_dir ${ARCS_SOURCE_DIR}/.. ABSOLUTE)
# addonのファイルを.cmakeから取得
include(${arcs_robot_dir}/equip/ARCS_EQUIP.cmake)

# CMakeLists.txtがあるフォルダをinclude
include_directories(${ARCS_SOURCE_DIR})

# -------実行ファイル生成部--------
# 実行ファイルARCSを生成
add_executable(
        ARCS
        ${ARCS_LIB_files}
        ${ARCS_SYS_files}
        ${ARCS_EQUIP_files}
        ${ARCS_CC_file} # ARCS.cc オンライン用
        ConstParams.hh
        ControlFunctions.cc
)

# ARCSに必要ライブラリをリンク
target_link_libraries(
        ARCS
        m
        ncursesw
        rt
        tinfo
        png
        z
        ${CMAKE_THREAD_LIBS_INIT}
)


# 実行ファイルARCS_offlineを生成
add_executable(
        ARCS_offline
        ${ARCS_LIB_files}
        ${ARCS_SYS_files}
        OfflineFunction.cc # ARCS.ccの代わりにオフライン
        ConstParams.hh
        ControlFunctions.cc
)

# ARCSに必要ライブラリをリンク
target_link_libraries(
        ARCS_offline
        m
        ncursesw
        rt
        tinfo
        png
        z
        ${CMAKE_THREAD_LIBS_INIT}
)