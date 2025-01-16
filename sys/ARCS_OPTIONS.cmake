# オプションに言語フィルタを適用する関数
function(add_language_specific_options var_name language options)
    set(language_specific_options "")
    foreach(option IN LISTS options)
        list(APPEND language_specific_options "$<$<COMPILE_LANGUAGE:${language}>:${option}>")
    endforeach()
    set(${var_name} ${language_specific_options} PARENT_SCOPE)
endfunction()

# 基本オプションの定義
set(CXX_OPTIONS_BASE
        -Wall
        -Weffc++
        -ftree-vectorize
        -march=native
        -fPIE
)

set(C_OPTIONS_BASE
        -Wall
        -ftree-vectorize
        -march=native
        -fPIE
)

# 言語固有のオプションを生成
add_language_specific_options(ARCS_COMMON_CXX_OPTIONS CXX "${CXX_OPTIONS_BASE}")
add_language_specific_options(ARCS_COMMON_C_OPTIONS C "${C_OPTIONS_BASE}")

# メインプロジェクト用の追加オプション
set(ARCS_MAIN_OPTIONS -DARCS_IN)

# ビルドタイプ毎のオプション
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -pipe")
set(CMAKE_C_FLAGS_RELEASE "-O2 -pipe")
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -Og")
set(CMAKE_C_FLAGS_DEBUG "-ggdb3 -Og")

# C++バージョン指定
set(CMAKE_CXX_STANDARD 17)