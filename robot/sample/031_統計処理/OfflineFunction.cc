//! @file OfflineFunction.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2024/07/23
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

// 追加のARCSライブラリをここに記述
#include "ArcsMatrix.hh"
#include "CsvManipulator.hh"
#include "Matrix.hh"
#include "Statistics.hh"

using namespace ARCS;
using namespace ArcsMatrix;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	// 行列全体の平均
	printf("\n◆ 行列全体の平均\n");
	constexpr ArcsMat<3,4> Am = {
		1,  1,  1,  1,
		2,  3, -2,  5,
		3, -1,  2,  4
	};
	disp(Am);
	printf("meanA = % g\n", Statistics::Mean(Am));
	constexpr double meanAx = Statistics::Mean(Am);	// コンパイル時に計算
	printf("meanAx = % g\n", meanAx);

	// 行列の横方向の平均
	printf("\n◆ 行列の横方向の平均\n");
	ArcsMat<3,1> vm;
	Statistics::MeanRow(Am, vm);
	disp(vm);
	constexpr auto vmx = Statistics::MeanRow(Am);	// コンパイル時に計算
	disp(vmx);

	// 行列の縦方向の平均
	printf("\n◆ 行列の縦方向の平均\n");
	ArcsMat<1,4> hm;
	Statistics::MeanColumn(Am, hm);
	disp(hm);
	constexpr auto hmx = Statistics::MeanColumn(Am);// コンパイル時に計算
	disp(hmx);

	// 行列全体の平均のテスト
	printf("◆ 行列全体の平均のテスト\n");
	Matrix<4,3> A = {
		1,  1,  1,  1,
		2,  3, -2,  5,
		3, -1,  2,  4
	};
	PrintMat(A);
	printf("Mean = % g\n\n", Statistics::Mean(A));
	
	// 行列の横方向の平均のテスト
	printf("◆ 行列の横方向の平均のテスト\n");
	Matrix<1,3> v;
	Statistics::MeanRow(A, v);
	PrintMat(v);
	
	// 行列の縦方向の平均のテスト
	printf("◆ 行列の縦方向の平均のテスト\n");
	Matrix<4,1> h;
	Statistics::MeanColumn(A, h);
	PrintMat(h);
	
	// 行列全体の分散のテスト
	printf("◆ 行列全体の分散のテスト\n");
	printf("Variance = % g\n\n", Statistics::Variance(A));
	
	// 行列の横方向の分散のテスト
	printf("◆ 行列の横方向の分散のテスト\n");
	Statistics::VarianceRow(A, v);
	PrintMat(v);
	
	// 行列の縦方向の分散のテスト
	printf("◆ 行列の縦方向の分散のテスト\n");
	Statistics::VarianceColumn(A, h);
	PrintMat(h);
	
	// 行列全体の標準偏差のテスト
	printf("◆ 行列全体の標準偏差のテスト\n");
	printf("Standard Deviation = % g\n\n", Statistics::StandardDeviation(A));
	
	// 行列の横方向の標準偏差のテスト
	printf("◆ 行列の横方向の標準偏差のテスト\n");
	Statistics::StandardDeviationRow(A, v);
	PrintMat(v);
	
	// 行列の縦方向の標準偏差のテスト
	printf("◆ 行列の縦方向の標準偏差のテスト\n");
	Statistics::StandardDeviationColumn(A, h);
	PrintMat(h);
	
	return EXIT_SUCCESS;	// 正常終了
}

