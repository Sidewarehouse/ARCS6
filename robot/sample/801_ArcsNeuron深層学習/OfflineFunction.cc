//! @file OfflineFunction.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2025/12/31
//! @author Yokokura, Yuki
//!
//! @par オフライン計算用のメインコード
//! - 「make offline」でコンパイルすると，いつものARCS制御用のコードは走らずに，
//!    このソースコードのみが走るようになる。
//! - ARCSライブラリはもちろんそのままいつも通り使用可能。
//! - 従って，オフラインで何か計算をしたいときに，このソースコードに記述すれば良い。
//!
// Copyright (C) 2011-2025 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// 基本のインクルードファイル
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <complex>
#include <array>
#include <any>

// 追加のARCSライブラリをここに記述
#include "ArcsMatrix.hh"
#include "ArcsNeuron.hh"
#include "CsvManipulator.hh"

using namespace ARCS;
using namespace ArcsMatrix;
using namespace ArcsNeuron;

// プロトタイプ宣言
void AutoDiffTestCode1(void);	//!< 自動微分テストコード1
void AutoDiffTestCode2(void);	//!< 自動微分テストコード2

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	AutoDiffTestCode1();	// 自動微分テストコード1
	//AutoDiffTestCode2();	// 自動微分テストコード2
	
	return EXIT_SUCCESS;	// 正常終了
}

//! @brief 自動微分テストコード1
void AutoDiffTestCode1(void){
	ArcsNeuStack gt;	// 自動微分スタック(勾配テープ)
	//ArcsNeu<double> x(&gt), W(&gt), V(&gt), b(&gt), y(&gt);	// エッジ変数
	ArcsNeu<double> x(&gt), y(&gt);
	ArcsNeu<float> b(&gt);

	// 異なる型同士での演算例
	double xx = 3.14, yy = 0;
	float  bb = 2.71;
	yy = xx + bb;
	printf("yy = xx + bb = %f + %f = %f\n", xx, bb, yy);

	// エッジ変数のメモリアドレス
	x.DispAddress("x");
	//W.DispAddress("W");
	//V.DispAddress("V");
	b.DispAddress("b");
	y.DispAddress("y");

	// エッジ変数に値をセット
	x = 3.14;
	//W = 10;
	//V = 5;
	b = 2.71;
	
	// 複合式の自動微分テスト
	y = x + b;	// 左辺値 + 左辺値
	//y = W*x;		// 左辺値*左辺値
	//y = b + W*x;	// 左辺値 + 右辺値(左辺値*左辺値)
	//y = W*x + b;	// 右辺値(左辺値*左辺値) + 左辺値
	//y = W*x + V*b;// 右辺値(左辺値*左辺値) + 右辺値(左辺値*左辺値)
	//y = W*x + W*b;// 右辺値(重複左辺値*左辺値) + 右辺値(重複左辺値*左辺値)
	//y = W*(x + b);// 左辺値*( 右辺値(左辺値 + 左辺値) )
	//y = (W + V)*x;// 右辺値(左辺値 + 左辺値)*左辺値 
	//y = (W + V)*(x + b);		// 右辺値(左辺値 + 左辺値)*右辺値(左辺値 + 左辺値)
	//y = W*x + W*b + V*x + V*b;	// 上記を展開した場合
	//y = x + b + x;		// 分岐される場合
	//y = W*x + x;			// 分岐される場合
	//y = ReLU(x);			// ReLU(左辺値)
	//y = ReLU(W*x + b);	// ReLU(右辺値)
	//y = V*ReLU(W*x + b);	// 左辺値*右辺値( ReLU(右辺値) )
	//y = ReLU(W*x + b) + x;	// スキップ接続がある場合
	
	// 自動微分スタックの表示
	gt.DispStack();			// 演算履歴の表示
	//gt.DispTempObjStack();	// 永続化された一時オブジェクト履歴の表示
	//gt.ClearGradient();		// 勾配をゼロ初期化
	y.SetGradient(1.23);		// 最終出力の勾配を設定(Loss)
	gt.UpdateGradient();	// 勾配を更新
	//gt.DispBackwardCalc();	// 逆方向計算の表示

	x.Disp("x");
	//W.Disp("W");
	//V.Disp("V");
	b.Disp("b");
	y.Disp("y");
	//gt.DispTempObjVar();	// 永続化された一時オブジェクトエッジ変数値の表示
}

//! @brief 自動微分テストコード2
void AutoDiffTestCode2(void){
	//ArcsNeuStack<std::any> gt;
	//ArcsNeu<int> x(&gt);
	/*
	std::any a;
	std::array<double, 3> b, y;
	b.fill(3.14);
	a = b;
	y = std::any_cast< std::array<double, 3> >(a);
	for(size_t i = 0; i < 3; ++i) printf("% f ", y[i]);
	*/
	
	// std::anyとポインタの組み合わせテスト
	double b = 3.14;
	double* c;
	c = &b;
	printf("*c = %f\n", *c);
	std::any a;
	a = c;
	printf("_a = %f\n", *std::any_cast<double*>(a));

	// std::anyとArcsMatポインタの組み合わせ
	ArcsMat<3,3> B(3.14);	// 行列の実体
	ArcsMat<3,3>* C;		// 行列のポインタ
	C = &B;					// ポインタにアドレスをリンク
	std::any A;				// どんな型でも入る変数A
	A = C;					// 行列のポインタ型をセット
	(*std::any_cast<ArcsMat<3,3>*>(A)).Disp();	// 行列のポインタ型として呼び出し
}
