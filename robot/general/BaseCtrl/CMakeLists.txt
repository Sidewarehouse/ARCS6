cmake_minimum_required(VERSION 3.16)
project(ARCS)

# pthreadを見つける
find_package(Threads REQUIRED)

# ARCSのルートディレクトリを取得
get_filename_component(arcs_root_dir ${ARCS_SOURCE_DIR}/../../.. ABSOLUTE)

# 共通オプションとメインプロジェクト用オプションを追加
include(${arcs_root_dir}/sys/ARCS_OPTIONS.cmake)
add_compile_options(
        ${ARCS_COMMON_CXX_OPTIONS}
        ${ARCS_COMMON_C_OPTIONS}
        ${ARCS_MAIN_OPTIONS}
)

# libを取得 (static リンク)
add_subdirectory(${arcs_root_dir}/lib lib)
# sysのファイルを.cmakeから取得
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
        ${ARCS_CC_file} # ARCS.cc オンライン用
        ${ARCS_SYS_files}
        ${ARCS_EQUIP_files}
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
        ARCS_LIB
)


# 実行ファイルARCS_offlineを生成
add_executable(
        ARCS_offline
        ${ARCS_SYS_files}
        ${ARCS_EQUIP_files}
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
        ARCS_LIB
)

# ソースコード行数カウント用のターゲット
add_custom_target(count_lines
        COMMAND wc -l
        ${ARCS_SOURCE_DIR}/*.cc
        ${ARCS_SOURCE_DIR}/*.hh
        ${arcs_robot_dir}/equip/*.cc
        ${arcs_robot_dir}/equip/*.hh
        ${arcs_robot_dir}/equip/c/*.c
        ${arcs_robot_dir}/equip/c/*.h
        ${arcs_root_dir}/lib/*.cc
        ${arcs_root_dir}/lib/*.hh
        ${arcs_root_dir}/sys/*.cc
        ${arcs_root_dir}/sys/*.hh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)