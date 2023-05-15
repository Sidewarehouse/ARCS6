//! @file OfflineCalculations.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2020/06/10
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
#include <chrono>

// 追加のARCSライブラリをここに記述
#include "Matrix.hh"
#include "CsvManipulator.hh"
#include "RandomGenerator.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	
	// 次数設定
	constexpr size_t N = 51;
	
	// 消費時間計測用
	double consmpt_time = 0;
	std::chrono::system_clock::time_point start_time, end_time;
	
	// 乱数行列の生成
	RandomGenerator Rnd(-1, 1);
	Matrix<N,N> A;
	Rnd.GetRandomMatrix(A);
	//PrintMat(A);
	
	// LU分解
	Matrix<N,N> L, U;
	Matrix<1,N,int> v;
	LU(A, L, U, v);
	//PrintMatrix(L, "% 6.3f");
	//PrintMatrix(U, "% 6.3f");
	
	// 普通の逆行列
	Matrix<N,N> Uinv;
	start_time = std::chrono::system_clock::now();
	Uinv = inv(U);
	end_time = std::chrono::system_clock::now();
	//PrintMatrix(Uinv, "% 6.3f");
	
	// 消費時間計算
	consmpt_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000.0);
	printf("普通の逆行列の消費時間 = %f [ms]\n\n", consmpt_time);
		
	// 上三角行列版の手抜き逆行列
	start_time = std::chrono::system_clock::now();
	inv_upper_tri(U, Uinv);
	end_time = std::chrono::system_clock::now();
	//PrintMatrix(Uinv, "% 6.3f");
	
	// 消費時間計算
	consmpt_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000.0);
	printf("上三角行列版の手抜き逆行列 = %f [ms]\n\n", consmpt_time);
	
	// 単位行列になるか一応チェック
	PrintMatrix(U*Uinv, "% 2.0f");
	
	return EXIT_SUCCESS;	// 正常終了
}

