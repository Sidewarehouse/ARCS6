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
	
	// 定義
	ArcsNeuStack<double> gt;	// 自動微分スタック(勾配テープ)
	ArcsNeu<double> x(&gt), W(&gt), V(&gt), b(&gt), y(&gt);	// エッジ変数
	
	// 変数のメモリアドレス
	x.DispAddress("x");
	W.DispAddress("W");
	V.DispAddress("V");
	b.DispAddress("b");
	y.DispAddress("y");

	// 変数値
	x = 3;
	W = 10;
	V = 5;
	b = 1.1;
	
	// 複合式の自動微分テスト
	//y = x + b;	// 左辺値 + 左辺値
	//y = W*x;		// 左辺値*左辺値
	//y = b + W*x;	// 左辺値 + 右辺値(左辺値*左辺値)
	y = W*x + b;	// 右辺値(左辺値*左辺値) + 左辺値
	//y = W*x + V*b;// 右辺値(左辺値*左辺値) + 右辺値(左辺値*左辺値)
	//y = W*x + W*b;// 右辺値(重複左辺値*左辺値) + 右辺値(重複左辺値*左辺値)
	//y = W*(x + b);// 左辺値*( 右辺値(左辺値 + 左辺値) )
	//y = (W + V)*x;// 右辺値(左辺値 + 左辺値)*左辺値 
	//y = (W + V)*(x + b);		// 右辺値(左辺値 + 左辺値)*右辺値(左辺値 + 左辺値)
	//y = W*x + W*b + V*x + V*b;	// 上記を展開した場合
	//y = x + b + x;	// 分岐される場合
	//y = W*x + x;	// スキップがある場合

	// 自動微分スタックの表示
	gt.DispStack();			// 演算履歴の表示
	gt.DispTempObjStack();	// 永続化された一時オブジェクトエッジ変数の表示
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

