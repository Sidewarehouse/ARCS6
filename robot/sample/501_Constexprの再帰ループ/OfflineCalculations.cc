//! @file OfflineCalculations.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2020/08/26
//! @author Yokokura, Yuki
//!
//! @par オフライン計算用のメインコード
//! - 「make offline」でコンパイルすると，いつものARCS制御用のコードは走らずに，
//!    このソースコードのみが走るようになる。
//! - ARCSライブラリはもちろんそのままいつも通り使用可能。
//! - 従って，オフラインで何か計算をしたいときに，このソースコードに記述すれば良い。
//!
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

// 基本のインクルードファイル
#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <array>
#include <complex>
#include <iostream>

// 追加のARCSライブラリをここに記述
#include "Matrix.hh"
#include "CsvManipulator.hh"

using namespace ARCS;

template <std::size_t ... Is>
auto iterate(std::index_sequence<Is...> const &){
	return std::make_tuple(Matrix<1,Is+1>{}...);
}

void PrintMatrices

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	
	auto matrices {
		iterate(std::make_index_sequence<3>{})
	};
	
	PrintMat(std::get<0>(matrices));
	PrintMat(std::get<1>(matrices));
	PrintMat(std::get<2>(matrices));
	
	return EXIT_SUCCESS;	// 正常終了
}

