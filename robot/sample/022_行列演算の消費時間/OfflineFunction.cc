//! @file OfflineFunction.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2024/07/22
//! @author Yokokura, Yuki
//!
//! @par オフライン計算用のメインコード
//! - 「make offline」でコンパイルすると，いつものARCS制御用のコードは走らずに，
//!    このソースコードのみが走るようになる。
//! - ARCSライブラリはもちろんそのままいつも通り使用可能。
//! - 従って，オフラインで何か計算をしたいときに，このソースコードに記述すれば良い。
//!
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// 基本のインクルードファイル
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <complex>
#include <array>
#include <chrono>

// 追加のARCSライブラリをここに記述
#include "ArcsMatrix.hh"
#include "CsvManipulator.hh"
#include "Matrix.hh"
#include "RandomGenerator.hh"

using namespace ARCS;
using namespace ArcsMatrix;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	constexpr size_t N = 12;	// 次数設定
	constexpr size_t K = 100;	// 測定する計算回数
	
	// 消費時間計測用
	double consmpt_time = 0;
	std::chrono::system_clock::time_point start_time, end_time;
	
	// 乱数行列の生成
	RandomGenerator Rnd(-1, 1);
	ArcsMat<N,N> A1;
	Rnd.GetRandomMatrix(A1);
	dispf(A1, "% 4.3f");
	Matrix<N,N> A2;
	Rnd.GetRandomMatrix(A2);
	PrintMatrix(A2, "% 4.3f");
	
	// ArcsMatrixでの計算
	ArcsMat<N,N> Y1;
	start_time = std::chrono::system_clock::now();
	for(size_t i = 1; i <= K; ++i) Y1 = expm(A1);
	end_time = std::chrono::system_clock::now();
	consmpt_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000.0);
	printf("ArcsMatrixでの消費時間 = %f [ms]\n\n", consmpt_time);
	
	// Matrixでの計算
	Matrix<N,N> Y2;
	start_time = std::chrono::system_clock::now();
	for(size_t i = 1; i <= K; ++i) Y2 = expm(A2, 13);
	end_time = std::chrono::system_clock::now();
	consmpt_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000.0);
	printf("Matrixでの消費時間     = %f [ms]\n\n", consmpt_time);

	return EXIT_SUCCESS;	// 正常終了
}
