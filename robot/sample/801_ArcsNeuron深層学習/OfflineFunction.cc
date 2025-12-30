//! @file OfflineFunction.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2024/06/25
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
#include "ArcsNeuron.hh"
#include "CsvManipulator.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	ArcsNeuStack<double> gt;
	ArcsNeu<double> x(&gt), W(&gt), V(&gt), b(&gt), y(&gt);
	
	x.DispAddress("x");
	W.DispAddress("W");
	V.DispAddress("V");
	b.DispAddress("b");
	y.DispAddress("y");
	
	x = 3;
	W = 10;
	V = 5;
	b = 1.1;
	
	//y = x + b;
	//y = x*b;
	//y = x + W*b;
	//y = W*x + b;
	//y = W*x + V*b;
	y = W*x + W*b;
	
	gt.DispStack();
	gt.DispTempObjStack();
	//*
	gt.ClearGradient();
	y.SetGradient(7);
	gt.UpdateGradient();
	gt.DispBackwardCalc();

	x.Disp("x");
	W.Disp("W");
	V.Disp("V");
	b.Disp("b");
	y.Disp("y");
	gt.DispTempObjStack();
	//*/

	return EXIT_SUCCESS;	// 正常終了
}

