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

// 追加のARCSライブラリをここに記述
#include "ArcsVarMatrix.hh"
//#include "CsvManipulator.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	
	// 行列宣言と値のセット
	printf("★★★★★★★ 行列宣言と値のセット\n");
	ArcsVarMat<double> A(3,3, {	// 宣言と同時に値をセットする場合
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	});
	ArcsVarMat B(3,3);			// 宣言したあとに、
	B.Set(						// 値をセットする場合
		5,  1, -3,
		3, -1,  7,
		4,  2, -5
	);
	auto C = A;					// 宣言と同時に既にある行列を代入する場合
	ArcsVarMat<int> Aint(2,2, {	// int型の行列を定義する場合
		 1, 3,
		-5, 7
	});
	//ArcsVarMat<std::string> Astr(2,2);	// 数値型以外だとエラー！ (アンコメントするとAssertion Failed)
	
	// 行列の表示
	printf("\n★★★★★★★ 行列の表示\n");
	A.Disp();				// 行列Aを表示
	B.Disp();				// 行列Bを表示
	C.Disp("% 6.4f");		// 表示の書式指定をして表示する場合
	Aint.Disp("% d");		// 整数行列Aintを表示
	A.DispAddress();		// 行列Aのメモリアドレスを表示

	// 行列サイズの変更 (ArcsMatでは実装不可)
	C.Resize(5,7);
	C.Disp();
	C.Resize(4,6);
	C.Disp();
	C.Resize(3,3);
	C.Disp();

	// 行列要素からの取得
	printf("\n★★★★★★★ 行列要素からの取得\n");
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	A.Get(
		a11, a12, a13,	// 各要素をそれぞれのdouble型変数に読み込む関数
		a21, a22, a23,
		a31, a32, a33
	);
	printf("a11 = % 3.1f, a12 = % 3.1f, a13 = % 3.1f\n", a11, a12, a13);
	printf("a21 = % 3.1f, a22 = % 3.1f, a23 = % 3.1f\n", a21, a22, a23);
	printf("a31 = % 3.1f, a32 = % 3.1f, a33 = % 3.1f\n", a31, a32, a33);
	
	// 行列の生配列
	printf("\n★★★★★★★ 行列の生配列\n");
	ArcsVarMat<double> J[3];
	for(size_t i = 0; i < 3; ++i){
		J[i].Resize(3,3);
		J[i].Set(
			i + 1, i + 1, i + 1,
			    2,     3,    -2,
			    3,    -1, i + 1
		);
		J[i].Disp();
	}
	
	// 行列のstd::array配列
	printf("\n★★★★★★★ 行列のstd::array配列\n");
	std::array<ArcsVarMat<double>, 3> K;
	for(size_t i = 0; i < 3; ++i){
		K.at(i).Resize(3,3);
		K.at(i).Set(
			i + 1, i + 1, i + 1,
			    2,     3,    -2,
			    3,    -1, i + 1
		);
		K.at(i).Disp();
	}
	
	// 定数行列 (コンパイル時定数行列はArcsMatでは実装可だったが、ArcsVarMatでは実装不可)
	printf("\n★★★★★★★ 定数行列 (コンパイル時定数行列は実装不可)\n");
	const ArcsVarMat Alpha(2,2, {	// 定数行列として定義
		1, 2,
		3, 4
	});
	Alpha.Disp("% d");

	// 縦ベクトルの[]オペレータによる要素アクセス
	printf("\n★★★★★★★ 縦ベクトルの[]オペレータによる要素アクセス\n");
	ArcsVarMat<double> alpha(5,1, {	// 縦ベクトルを定義
		1,
		2,
		3,
		4,
		5
	});
	alpha.Disp();
	printf("alpha[3] = % g\n\n", alpha[3]);	// 縦ベクトルの3番目を表示
	alpha[3] = 6;			// 縦ベクトルの3番目に 6 を書き込む
	alpha.Disp();			// 縦ベクトルの3番目が変わっている
	
	// 行列の()オペレータによる要素アクセス(サイズチェック無し版)
	printf("\n★★★★★★★ 行列の()オペレータによる要素アクセス(サイズチェック無し版)\n");
	printf("B(1,1) = % g\n", B(1,1));	// 引数の順番は、B(縦方向, 横方向) であり、
	printf("B(1,2) = % g\n", B(1,2));	// MATLABの定義と同一である
	printf("B(1,3) = % g\n", B(1,3));
	printf("B(2,1) = % g\n", B(2,1));
	printf("B(2,2) = % g\n", B(2,2));
	printf("B(2,3) = % g\n", B(2,3));
	printf("B(3,1) = % g\n", B(3,1));
	printf("B(3,2) = % g\n", B(3,2));
	printf("B(3,3) = % g\n", B(3,3));
	B(2,3) = 0;	// 2行3列目に 0 を書き込む
	B.Disp();	// 変わっている
	B(2,3) = 7;	// 2行3列目に 7 を書き込む
	B.Disp();	// 元に戻っている
	
	// 行列の()オペレータによる要素アクセス(サイズチェック可能版)
	printf("\n★★★★★★★ 行列の()オペレータによる要素アクセス(サイズチェック可能版)\n");
	printf("B(2,3,true) = % g\n", B(2,3,true));		// 何も問題ない場合 (trueを付けるとサイズチェック有り)
	//printf("B(2,4,true) = % g\n", B(2,4,true));	// ←サイズエラー！
	//printf("B(2,4,true) = % g\n", B(2,4,false));	// ←運が悪いとSIGSEGV (falseにするとサイズチェック無し)
	B(2,3,true) = 0;	// 2行3列目に 0 を書き込む
	B.Disp();			// 変わっている
	B(2,3,true) = 7;	// 2行3列目に 7 を書き込む
	B.Disp();			// 元に戻っている
	//B(2,4,true) = 7;	// ←サイズエラー！
	//B(2,4,false) = 7;	// ←運が悪いとSIGSEGV
	
	// 代入演算子
	printf("\n★★★★★★★ 代入演算子\n");
	C = B;
	C.Disp();
	//C = Alpha;	// ←サイズエラー！
/*	
	// 演算子
	printf("\n★★★★★★★ 演算子\n");
	C = +A;					// 単項加算
	printf("C = +A \n");
	disp(C);
	C = -A;					// 単項減算
	printf("C = -A\n");
	disp(C);
	C = A + B;				// 行列加算
	printf("C = A + B\n");
	disp(C);
	C = A + 3;				// 行列スカラー加算
	printf("C = A + 3\n");
	disp(C);
	C = A - 3;				// 行列スカラー減算
	printf("C = A - 3\n");
	disp(C);
	C = A - B;				// 行列減算
	printf("C = A - B\n");
	disp(C);
	C = A*B;				// 行列乗算
	printf("C = A*B\n");
	disp(C);
	C = A*2;				// 行列スカラー乗算
	printf("C = A*2\n");
	disp(C);
	C = A/2;				// 行列スカラー除算
	printf("C = A/2\n");
	disp(C);
	C += A;					// 行列加算代入
	printf("C += A\n");
	disp(C);
	C += 2;					// 行列スカラー加算代入
	printf("C += 2\n");
	disp(C);
	C -= A;					// 行列減算代入
	printf("C -= A\n");
	disp(C);
	C -= 1;					// 行列スカラー減算代入
	printf("C -= 1\n");
	disp(C);
	C *= A;					// 正方行列乗算代入
	printf("C *= A\n");
	disp(C);
	C *= 2;					// 行列スカラー乗算代入
	printf("C *= 2\n");
	disp(C);
	C /= 2;					// 行列スカラー除算代入
	printf("C /= 2\n");
	disp(C);
	C = A^0;				// 正方行列の0のべき乗
	printf("C = A^0\n");
	disp(C);
	C = A^1;				// 正方行列の1のべき乗
	printf("C = A^1\n");
	disp(C);
	C = A^2;				// 正方行列の2のべき乗
	printf("C = A^2\n");
	disp(C);
	C = A^3;				// 正方行列の3のべき乗
	printf("C = A^3\n");
	disp(C);
	C = A & B;				// アダマール積
	printf("C = A & B\n");
	disp(C);
	C = A % B;				// 要素ごとの除算
	printf("C = A %% B\n");
	disp(C);
	C = 3 + A;				// スカラー行列加算
	printf("C = 3 + A\n");
	disp(C);
	C = 3 - A;				// スカラー行列減算
	printf("C = 3 - A\n");
	disp(C);
	C = 2*A;				// スカラー行列乗算
	printf("C = 2*A\n");
	disp(C);
*/
	return EXIT_SUCCESS;	// 正常終了
}

