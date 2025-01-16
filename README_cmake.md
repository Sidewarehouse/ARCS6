# ARCS-PUBLIC

## CMakeでビルド

[Makefileの時と基本的な考え方は同一](https://www.sidewarehouse.net/arcs6/compile.html)

general-BaseCtrlの場合

```bash
cd path/to/ARCS6
cd robot/general/BaseCtrl
mkidr build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j ARCS # or ARCS_offline
```

で実行ファイルARCSがコンパイルされる．

### オプション

```bash
cmake -DCMAKE_BUILD_TYPE=Debug .. #GDBデバッグ用にビルド
cmake -DCMAKE_BUILD_TYPE=ASM .. #アセンブリリストを出力
make VERBOSE=1 #で実際に実行されているコマンドが見える
```

### CLionでの使い方

JetBrains CLionを用いてUbuntu Serverがインストールされている制御用PCにフルリモートモードで開発するために，CMake対応にした．

通常の設定方法だとBaseCtrl以下のファイルしか同期されないため，CLionのPreferences-Build, Execution, Deployment-MappingsのLocal path:をARCS6に変更する．

ロボットや制御を変更するときは，CMakeLists.txtを右クリックしてLoad CMake Projectする．

## ライブラリとして使用
このように記述することで，ARCS6をライブラリとして使用できる．

```cmake
include(FetchContent)
FetchContent_Declare(
        ARCS
        GIT_REPOSITORY https://github.com/Sidewarehouse/ARCS6.git
        GIT_TAG main
)
FetchContent_MakeAvailable(ARCS)
add_subdirectory(${arcs_SOURCE_DIR}/lib)
add_executable(BIN_NAME
        main.cpp
)
target_link_libraries(BIN_NAME
        PRIVATE
        ARCS_LIB
)
```
