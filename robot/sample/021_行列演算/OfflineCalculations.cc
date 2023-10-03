//! @file OfflineCalculations.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2022/08/05
//! @author Yokokura, Yuki
//!
//! @par オフライン計算用のメインコード
//! - 「make offline」でコンパイルすると，いつものARCS制御用のコードは走らずに，
//!    このソースコードのみが走るようになる。
//! - ARCSライブラリはもちろんそのままいつも通り使用可能。
//! - 従って，オフラインで何か計算をしたいときに，このソースコードに記述すれば良い。
//!
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

// 基本のインクルードファイル
#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <array>
#include <complex>
#include <cmath>

// 追加のARCSライブラリをここに記述
#include "ArcsMatrix.hh"
//#include "CsvManipulator.hh"

using namespace ARCS;
using namespace ArcsMatrix;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	
	// 行列宣言と値のセット
	printf("★★★★★★★ 行列宣言と値のセット\n");
	ArcsMat<3,3> A = {		// 宣言と同時に値をセットする場合
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	};
	ArcsMat<3,3> B;			// 宣言したあとに、
	B.Set(					// 値をセットする場合
		5,  1, -3,
		3, -1,  7,
		4,  2, -5
	);
	auto C = A;				// 宣言と同時に既にある行列を代入する場合
	ArcsMat<2,3> Pi(M_PI);	// 宣言と同時に設定値で全て埋める場合
	
	// 行列の表示
	printf("\n★★★★★★★ 行列の表示\n");
	dispsize(A);				// 行列Aのサイズを表示
	dispmat(A);					// 行列Aを表示
	dispmat(B);					// 行列Bを表示
	dispmatfmt(C, "% 6.4f");	// 表示の書式指定をして表示する場合
	dispmatfmt(Pi, "% 16.15f");	// 行列Piを書式指定して表示
	A.DispAddress();			// 行列Aのメモリアドレスを表示
	
	// 行列のサイズの取得
	printf("H = %zu, W = %zu\n", A.GetHeight(), A.GetWidth());
	
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
	ArcsMat<3,3> J[3];
	for(size_t i = 0; i < 3; ++i){
		J[i].Set(
			i + 1, i + 1, i + 1,
			    2,     3,    -2,
			    3,    -1, i + 1
		);
		dispmat(J[i]);
	}
	
	// 行列のstd::array配列
	printf("\n★★★★★★★ 行列のstd::array配列\n");
	std::array<ArcsMat<3,3>, 3> K;
	for(size_t i = 0; i < 3; ++i){
		K.at(i).Set(
			i + 1, i + 1, i + 1,
			    2,     3,    -2,
			    3,    -1, i + 1
		);
		dispmat(K.at(i));
	}
	
	// 定数行列，コンパイル時定数行列と初期化
	printf("\n★★★★★★★ 定数行列，コンパイル時定数行列と初期化\n");
	const ArcsMat<2,2> Alpha = {	// 定数行列として定義
		1, 2,
		3, 4
	};
	constexpr ArcsMat<3,3> Ax = {	// コンパイル時定数行列として定義
		5, 6, 7,
		8, 9, 0,
		1, 2, 3
	};
	dispmat(Alpha);
	dispmat(Ax);
	
	// 縦ベクトルの[]オペレータによる要素アクセス
	printf("\n★★★★★★★ 縦ベクトルの[]オペレータによる要素アクセス\n");
	ArcsMat<5,1> alpha = {	// 縦ベクトルを定義
		1,
		2,
		3,
		4,
		5
	};
	dispsize(alpha);
	dispmat(alpha);
	printf("alpha[3] = % g\n", alpha[3]);	// 縦ベクトルの3番目を表示
	alpha[3] = 6;			// 縦ベクトルの3番目に 6 を書き込む
	dispmat(alpha);			// 縦ベクトルの3番目が変わっている
	
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
	dispmat(B);	// 変わっている
	B(2,3) = 7;	// 2行3列目に 7 を書き込む
	dispmat(B);	// 元に戻っている
	
	// 行列の()オペレータによる要素アクセス(サイズチェック可能版)
	printf("\n★★★★★★★ 行列の()オペレータによる要素アクセス(サイズチェック可能版)\n");
	printf("B(2,3,true) = % g\n", B(2,3,true));		// 何も問題ない場合 (trueを付けるとサイズチェック有り)
	//printf("B(2,4,true) = % g\n", B(2,4,true));	// ←サイズエラー！
	//printf("B(2,4,true) = % g\n", B(2,4,false));	// ←運が悪いとSIGSEGV (falseにするとサイズチェック無し)
	B(2,3,true) = 0;	// 2行3列目に 0 を書き込む
	dispmat(B);			// 変わっている
	B(2,3,true) = 7;	// 2行3列目に 7 を書き込む
	dispmat(B);			// 元に戻っている
	//B(2,4,true) = 7;	// ←サイズエラー！
	//B(2,4,false) = 7;	// ←運が悪いとSIGSEGV
	
	// 代入演算子
	printf("\n★★★★★★★ 代入演算子\n");
	C = B;
	dispmat(C);
	//C = Alpha;	// ←サイズエラー！
	
	// 演算子
	printf("\n★★★★★★★ 演算子\n");
	C = +A;					// 単項加算
	printf("C = +A \n");
	dispmat(C);
	C = -A;					// 単項減算
	printf("C = -A\n");
	dispmat(C);
	C = A + B;				// 行列加算
	printf("C = A + B\n");
	dispmat(C);
	C = A + 3;				// 行列スカラー加算
	printf("C = A + 3\n");
	dispmat(C);
	C = A - 3;				// 行列スカラー減算
	printf("C = A - 3\n");
	dispmat(C);
	C = A - B;				// 行列減算
	printf("C = A - B\n");
	dispmat(C);
	C = A*B;				// 行列乗算
	printf("C = A*B\n");
	dispmat(C);
	C = A*2;				// 行列スカラー乗算
	printf("C = A*2\n");
	dispmat(C);
	C = A/2;				// 行列スカラー除算
	printf("C = A/2\n");
	dispmat(C);
	C += A;					// 行列加算代入
	printf("C += A\n");
	dispmat(C);
	C += 2;					// 行列スカラー加算代入
	printf("C += 2\n");
	dispmat(C);
	C -= A;					// 行列減算代入
	printf("C -= A\n");
	dispmat(C);
	C -= 1;					// 行列スカラー減算代入
	printf("C -= 1\n");
	dispmat(C);
	C *= A;					// 正方行列乗算代入
	printf("C *= A\n");
	dispmat(C);
	C *= 2;					// 行列スカラー乗算代入
	printf("C *= 2\n");
	dispmat(C);
	C = A^0;				// 正方行列の0のべき乗
	printf("C = A^0\n");
	dispmat(C);
	C = A^1;				// 正方行列の1のべき乗
	printf("C = A^1\n");
	dispmat(C);
	C = A^2;				// 正方行列の2のべき乗
	printf("C = A^2\n");
	dispmat(C);
	C = A^3;				// 正方行列の3のべき乗
	printf("C = A^3\n");
	dispmat(C);
	C = A & B;				// アダマール積
	printf("C = A & B\n");
	dispmat(C);
	C = A % B;				// アダマール除算
	printf("C = A %% B\n");
	dispmat(C);
	C = 3 + A;				// スカラー行列加算
	printf("C = 3 + A\n");
	dispmat(C);
	C = 3 - A;				// スカラー行列減算
	printf("C = 3 - A\n");
	dispmat(C);
	C = 2*A;				// スカラー行列乗算
	printf("C = 2*A\n");
	dispmat(C);
	
	// 下記はサイズエラー等の場合(アンコメントするとAssertion Failed)
	//C = A + Alpha;	// ←サイズエラー！
	//C = A - Alpha;	// ←サイズエラー！
	//C = A*Alpha;		// ←サイズエラー！
	//C += Alpha;		// ←サイズエラー！
	//C -= Alpha;		// ←サイズエラー！
	//C *= Alpha;		// ←サイズエラー！
	//alpha *= alpha;	// ←正方行列エラー！
	//alpha = alpha^2;	// ←正方行列エラー！
	//C = A & Alpha;	// ←サイズエラー！
	//C = A % Alpha;	// ←サイズエラー！
	//C = 1/A;			// ←この演算子の使い方は使用禁止！
	
	// すべての要素を埋める
	printf("\n★★★★★★★ すべての要素を埋める\n");
	C.FillAll(3.9);	// すべての要素を3.9で埋める
	dispmat(C);
	C.FillAllZero();// すべての要素をゼロで埋める
	dispmat(C);
	C = B;
	
	// 1次元std::array配列との相互変換
	printf("\n★★★★★★★ 1次元std::array配列との相互変換\n");
	std::array<double, 10> v10array = {100,200,300,400,500,600,700,800,900,1000};
	ArcsMat<10,1> v10;
	v10.LoadArray(v10array);	// std::arrayから読み込む
	dispmat(v10);
	v10[10] = 3939;
	v10.StoreArray(v10array);	// std::arrayに書き込む
	printf("v10array = ");
	for(size_t i = 0; i < 10; ++i) printf("% 5.0f ", v10array.at(i));
	printf("\n");
	//A.LoadArray(v10array);	// ←サイズエラー！(アンコメントするとAssertion Failed)
	//std::array<std::string, 10> v10chars = {"AAA"};
	//v10.LoadArray(v10chars);	// ←型変換エラー！(アンコメントするとAssertion Failed)
	
	// 2次元std::array配列との相互変換
	printf("\n★★★★★★★ 2次元std::array配列との相互変換\n");
	std::array<std::array<double, 3>, 2> D32array = {{
		{1, 2, 3},
		{4, 5, 6}
	}};
	ArcsMat<3,2> D;
	D.LoadArray(D32array);		// std::arrayから読み込む(ただし転置した状態となる)
	dispmat(D);
	D(3,2) = 39;
	D.StoreArray(D32array);		// std::arrayに書き込む(ただし転置した状態となる)
	printf("D32array:1 = ");
	for(size_t i = 0; i < 3; ++i) printf("% 5.0f ", D32array[0][i]);
	printf("\nD32array:2 = ");
	for(size_t i = 0; i < 3; ++i) printf("% 5.0f ", D32array[1][i]);
	printf("\n");
	//std::array<std::array<double, 4>, 3> D43array;
	//D.LoadArray(D43array);	// ←サイズエラー！(アンコメントするとAssertion Failed)
	
	// ベクトルの抜き出し
	printf("\n★★★★★★★ ベクトルの抜き出し\n");
	dispmat(C);
	ArcsMat<2,1> v2;
	C.GetVerticalVec(v2,2,3);			// 2行3列目を先頭として高さ2の縦ベクトルを抜き出す (引数渡し版)
	dispmat(v2);
	dispmat(C.GetVerticalVec<2>(2,3));	// 2行3列目を先頭として高さ2の縦ベクトルを抜き出す (戻り値渡し版)
	ArcsMat<1,2> w2;
	C.GetHorizontalVec(w2,2,2);			// 2行2列目を先頭として幅2の横ベクトルを抜き出す (引数渡し版)
	dispmat(w2);	
	dispmat(C.GetHorizontalVec<2>(2,2));	// 2行2列目を先頭として幅2の横ベクトルを抜き出す (戻り値渡し版)
	//dispmat(C.GetVerticalVec<2>(3,3));	// ←はみ出るエラー！(アンコメントするとAssertion Failed)
	//dispmat(C.GetHorizontalVec<2>(2,3));	// ←はみ出るエラー！(アンコメントするとAssertion Failed)
	constexpr auto v2x = Ax.GetVerticalVec<2>(1,3);		// コンパイル時に抜き出す
	constexpr auto w2x = Ax.GetHorizontalVec<2>(3,1);	// コンパイル時に抜き出す
	dispmat(v2x);
	dispmat(w2x);
	
	// ベクトルの埋め込み
	printf("\n★★★★★★★ ベクトルの埋め込み\n");
	v2.Set(
		39,
		93
	);
	C.SetVerticalVec(v2,2,3);		// 2行3列目を先頭として縦ベクトルを埋め込む
	dispmat(C);
	w2.Set(39, 93);
	C.SetHorizontalVec(w2,2,1);		// 2行1列目を先頭として横ベクトルを埋め込む
	dispmat(C);
	//C.SetVerticalVec(v2,3,1);		// ←はみ出るエラー！(アンコメントするとAssertion Failed)
	//C.SetHorizontalVec(w2,1,3);	// ←はみ出るエラー！(アンコメントするとAssertion Failed)

	// 零行列と1行列と単位行列
	printf("\n★★★★★★★ 零行列と1行列と単位行列\n");
	auto O = zeros<5,5>();	// 零行列
	auto l = ones<5,5>();	// 1行列
	auto I = eye<5,5>();	// 単位行列
	dispmat(O);
	dispmat(l);
	dispmat(I);
	constexpr auto Ox = zeros<3,5>();	// コンパイル時に零行列を生成
	constexpr auto lx = ones<3,5>();	// コンパイル時に1行列を生成
	constexpr auto Ix = eye<3,3>();		// コンパイル時に単位行列を生成
	dispmat(Ox);
	dispmat(lx);
	dispmat(Ix);
	
	// 単調増加ベクトル
	printf("\n★★★★★★★ 単調増加ベクトル\n");
	auto r = ramp<5,1>();				// 単調増加ベクトル
	constexpr auto rx = ramp<5,1>();	// コンパイル時に単調増加ベクトルを生成
	dispmat(r);
	dispmat(rx);
	
	// 転置行列
	printf("\n★★★★★★★ 転置行列\n");
	ArcsMat<2,3> Dt;
	dispmat(D);
	tp(D, Dt);		// Dの転置 (引数渡し版)
	dispmat(Dt);
	dispmat(tp(D));	// Dの転置 (戻り値渡し版)
	dispmat(Ax);
	constexpr auto Atx = tp(Ax);		// コンパイル時に転置行列を生成
	dispmat(Atx);
	//tp(A, Dt);	// ←サイズエラー！(アンコメントするとAssertion Failed)
	//A = tp(Dt);	// ←サイズエラー！(アンコメントするとAssertion Failed)
	
	// 列操作関連の関数
	printf("\n★★★★★★★ 列操作関連の関数\n");
	ArcsMat<3,1> vd;
	getcolumn(D, vd, 2);			// 行列Dの2列目を縦ベクトルとして抽出 (引数渡し版)
	dispmat(vd);
	dispmat(getcolumn(D, 1));		// 行列Dの1列目を縦ベクトルとして抽出 (戻り値渡し版)
	constexpr auto vax = getcolumn(Ax, 2);				// コンパイル時に縦ベクトルを抽出
	dispmat(vax);
	setcolumn(D, vax, 1);			// 行列Dの1列目に縦ベクトルを埋め込む (引数渡し版)
	dispmat(D);
	dispmat(setcolumn(vd, 2, D));	// 行列Dの2列目に縦ベクトルを埋め込む (戻り値渡し版)
	constexpr ArcsMat<3,4> Ex = {
		0, 1, 2, 3,
		4, 5, 6, 7,
		8, 9, 0, 1
	};
	constexpr ArcsMat<3,1> vex = {
		39,
		39,
		24
	};
	dispmat(Ex);
	dispmat(vex);
	constexpr auto Ex2 = setcolumn(vex, 3, Ex);			// コンパイル時に縦ベクトルを埋め込む
	dispmat(Ex2);
	ArcsMat<3,4> E = Ex;
	dispmat(E);
	swapcolumn(E, 1, 3);				// 行列Eの1列目と3列目を入れ替える (引数渡し版)
	dispmat(E);
	dispmat(swapcolumn(1, 3, E));		// 行列Eの1列目と3列目を入れ替える (戻り値渡し版)
	constexpr auto Ex3 = swapcolumn(1, 3, Ex2);			// コンパイル時に列を入れ替える
	dispmat(Ex3);
	fillcolumn(E, 39, 4, 2, 3);			// 行列Eの4列目2行目～3行目を39で埋める (引数渡し版)
	dispmat(E);
	dispmat(fillcolumn(1, 4, 2, 3, E));	// 行列Eの4列目2行目～3行目を1で埋める (戻り値渡し版)
	constexpr auto Ex4 = fillcolumn(39, 3, 1, 2, Ex3);	// コンパイル時に指定列を値で埋める
	dispmat(Ex4);
	ArcsMat<1,4,size_t> wor = {2, 1, 4, 3};
	dispmat(wor);
	E = ordercolumn(E, wor);		// 行列Eをworが昇順になるように列を並び替える (戻り値渡し版のみ)
	dispmat(E);
	E = ordercolumn(E, wor);		// 2回やると元に戻る
	dispmat(E);
	ordercolumn_and_vec(E, wor);	// 行列Eをworが昇順になるようにEとworの両方の列を並び替える (引数渡し版のみ)
	dispmat(E);
	dispmat(wor);
	constexpr ArcsMat<1,4,size_t> worx = {2, 1, 4, 3};
	constexpr auto Ex5 = ordercolumn(Ex, worx);			// コンパイル時に並び替える
	dispmat(Ex5);
	ArcsMat<1,4> e5;	
	sumcolumn(Ex5, e5);				// 各列の総和を計算する (引数渡し版)
	dispmat(e5);
	dispmat(sumcolumn(Ex5));		// 各列の総和を計算する (戻り値渡し版)
	constexpr auto ex5 = sumcolumn(Ex5);				// コンパイル時に総和を計算する
	dispmat(ex5);
	
	// 行操作関連の関数
	printf("\n★★★★★★★ 行操作関連の関数\n");
	ArcsMat<1,2> wd;
	getrow(D, wd, 2);				// 行列Dの2行目を横ベクトルとして抽出 (引数渡し版)
	dispmat(wd);
	dispmat(getrow(D, 3));			// 行列Dの3行目を横ベクトルとして抽出 (戻り値渡し版)
	constexpr auto wax = getrow(Ax, 3);					// コンパイル時に横ベクトルを抽出
	dispmat(wax);
	setrow(D, wd, 1);				// 行列Dの1行目に横ベクトルを埋め込む (引数渡し版)
	dispmat(D);
	dispmat(setrow(wd, 2, D));		// 行列Dの2行目に横ベクトルを埋め込む (戻り値渡し版)
	constexpr ArcsMat<4,3> Etx = {
		0, 1, 2,
		3, 4, 5,
		6, 7, 8,
		9, 0, 1
	};
	constexpr ArcsMat<1,3> wetx = {39, 39, 24};
	dispmat(Etx);
	dispmat(wetx);
	constexpr auto Etx2 = setrow(wetx, 3, Etx);			// コンパイル時に横ベクトルを埋め込む
	dispmat(Etx2);
	ArcsMat<4,3> Et = Etx;
	dispmat(Et);
	swaprow(Et, 1, 3);					// 行列Etの1行目と3行目を入れ替える (引数渡し版)
	dispmat(Et);
	dispmat(swaprow(1, 3, Et));			// 行列Etの1行目と3行目を入れ替える (戻り値渡し版)
	constexpr auto Etx3 = swaprow(1, 3, Etx2);			// コンパイル時に行を入れ替える
	dispmat(Etx3);
	fillrow(Et, 39, 4, 2, 3);			// 行列Etの4行目2列目～3列目を39で埋める (引数渡し版)
	dispmat(Et);
	dispmat(fillrow(1, 4, 2, 3, Et));	// 行列Etの4行目2列目～3列目を1で埋める (戻り値渡し版)
	constexpr auto Etx4 = fillrow(39, 3, 1, 2, Etx3);	// コンパイル時に指定行を値で埋める
	dispmat(Etx4);
	ArcsMat<4,1,size_t> vor = {
		2,
		1,
		4,
		3
	};
	dispmat(vor);
	dispmat(orderrow(Et, vor));			// 行列Etをvorが昇順になるように行を並び替える (戻り値渡し版のみ)
	dispmat(orderrow(Et, vor));			// 2回やると元に戻る
	orderrow_and_vec(Et, vor);			// 行列Etをvorが昇順になるようにEtとvorの両方の行を並び替える (引数渡し版のみ)
	dispmat(Et);
	dispmat(vor);
	constexpr ArcsMat<4,1,size_t> vorx = {
		2,
		1,
		4,
		3
	};
	constexpr auto Etx5 = orderrow(Etx, vorx);			// コンパイル時に並び替える
	dispmat(Etx5);
	ArcsMat<4,1> et5;
	sumrow(Etx5, et5);				// 各々の行の総和を計算する (引数渡し版)
	dispmat(et5);
	dispmat(sumrow(Etx5));			// 各々の行の総和を計算する (戻り値渡し版)
	constexpr auto etx5 = sumrow(Etx5);					// コンパイル時に総和を計算する
	dispmat(etx5);
	
	// 小行列操作関連の関数
	printf("\n★★★★★★★ 小行列操作関連の関数\n");
	constexpr ArcsMat<7,6> Fx = {
		10, 20, 30, 40, 50, 60,
		70, 80, 90, 10, 11, 12,
		13, 14, 15, 16, 17, 18,
		19, 20, 21, 22, 23, 24,
		25, 26, 27, 28, 29, 30,
		31, 32, 33, 34, 35, 36,
		37, 38, 39, 40, 41, 42
	};
	auto F = Fx;
	dispmat(F);
	constexpr ArcsMat<4,1> vfx = {
		99,
		88,
		77,
		66
	};
	auto vf = vfx;
	ArcsMat<3,1> vfg;
	getvvector(F, vfg, 3, 2);			// 行列Fの先頭位置(3,2)から縦ベクトルを抽出する (引数渡し版)
	dispmat(vfg);
	dispmat(getvvector<3>(F, 2, 3));	// 行列Fの先頭位置(2,3)から縦ベクトルを抽出する (戻り値渡し版)
	constexpr auto vfgx = getvvector<5>(Fx, 2, 5);		// コンパイル時に縦ベクトルを抽出
	dispmat(vfgx);
	setvvector(F, vf, 3, 5);			// 行列Fの先頭位置(3,5)に縦ベクトルvfを上書きする (引数渡し版)
	dispmat(F);
	dispmat(setvvector(vf, 3, 4, F));	// 行列Fの先頭位置(3,4)に縦ベクトルvfを上書きする (戻り値渡し版)
	constexpr auto Fx2 = setvvector(vfx, 3, 4, Fx);		// コンパイル時に縦ベクトルを上書き
	dispmat(Fx2);
	constexpr ArcsMat<1,4> wfx = { 55, 44, 33, 22 };
	auto wf = wfx;
	ArcsMat<1,3> wfg;
	gethvector(F, wfg, 2, 1);			// 行列Fの先頭位置(2,1)から横ベクトルを抽出する (引数渡し版)
	dispmat(wfg);
	dispmat(gethvector<3>(F, 5, 2));	// 行列Fの先頭位置(5,2)から3列の横ベクトルを抽出する (戻り値渡し版)
	constexpr auto wfgx = gethvector<5>(Fx, 6, 2);		// コンパイル時に横ベクトルを抽出
	dispmat(wfgx);
	sethvector(F, wf, 2, 3);			// 行列Fの先頭位置(2,3)に横ベクトルwfを上書きする (引数渡し版)
	dispmat(F);
	dispmat(sethvector(wf, 1, 3, F));	// 行列Fの先頭位置(1,3)に横ベクトルwfを上書きする (戻り値渡し版)
	constexpr auto Fx4 = sethvector(wfx, 2, 3, Fx2);	// コンパイル時に横ベクトルを上書き
	dispmat(Fx4);
	ArcsMat<4,3> F43;
	getsubmatrix(F, F43, 2, 1);			// 行列Fの左上位置(2,1)から4×3の小行列を抽出する (引数渡し版)
	dispmat(F43);
	dispmat(( getsubmatrix<4,3>(F, 2, 2) ));	// 行列Fの左上位置(2,2)から4×3の小行列を抽出する (戻り値渡し版)
	constexpr auto Fx5 = getsubmatrix<2,3>(Fx4, 3, 3);	// コンパイル時に小行列を抽出する
	dispmat(Fx5);
	F43.Set(
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		0, 1, 2
	);
	setsubmatrix(F, F43, 2, 4);			// 行列Fの左上位置(2,4)に4×3の行列を書き込む (引数渡し版)
	dispmat(F);
	dispmat(setsubmatrix(F43, 2, 1, F));// 行列Fの左上位置(2,1)に4×3の行列を書き込む (戻り値渡し版)
	constexpr ArcsMat<2,3> F23x = {
		100, 200, 300,
		400, 500, 600
	};
	constexpr auto Fx6 = setsubmatrix(F23x, 3, 3, Fx4);	// コンパイル時に小行列を書き込む
	dispmat(Fx6);
	
	// シフト関連の関数
	printf("\n★★★★★★★ シフト関連の関数\n");
	F = Fx;
	auto G = ones<7,6>();
	dispmat(F);
	shiftup(F, G);			// 行列Fを1行分上にシフトして行列Gとして出力 (引数渡し版)
	dispmat(G);
	shiftup(G, G, 3);		// 行列Gを3行分上にシフト (引数渡し版)
	dispmat(G);
	dispmat(shiftup(F));	// 行列Fを1行分上にシフトして出力 (戻り値渡し版)
	dispmat(shiftup(F, 4));	// 行列Fを4行分上にシフトして出力 (戻り値渡し版)
	constexpr auto Gx1 = shiftup(Fx, 4);				// コンパイル時に上にシフト
	dispmat(Gx1);
	dispmat(F);
	shiftdown(F, G);		// 行列Fを1行分下にシフトして行列Gとして出力 (引数渡し版)
	dispmat(G);
	shiftdown(F, G, 3);		// 行列Fを3行分下にシフトして行列Gとして出力 (引数渡し版)
	dispmat(G);
	F = Fx;
	dispmat(F);
	dispmat(shiftdown(F));		// 行列Fを1行分下にシフトして出力 (戻り値渡し版)
	dispmat(shiftdown(F, 4));	// 行列Fを4行分下にシフトして出力 (戻り値渡し版)
	constexpr auto Gx2 = shiftdown(Fx, 2);				// コンパイル時に下にシフト
	dispmat(Gx2);
	F = Fx;
	shiftleft(F, G);			// 行列Fを1列分左にシフトして行列Gとして出力 (引数渡し版)
	dispmat(G);
	shiftleft(F, G, 3);			// 行列Fを3列分左にシフトして行列Gとして出力 (引数渡し版)
	dispmat(G);
	dispmat(shiftleft(F));		// 行列Fを1列分左にシフトして出力 (戻り値渡し版)
	dispmat(shiftleft(F, 2));	// 行列Fを2列分左にシフトして出力 (戻り値渡し版)
	constexpr auto Gx3 = shiftleft(Fx, 4);				// コンパイル時に左にシフト
	dispmat(Gx3);
	F = Fx;
	shiftright(F, G);			// 行列Fを1列分右にシフトして行列Gとして出力 (引数渡し版)
	dispmat(G);
	shiftright(F, G, 3);		// 行列Fを3列分右にシフトして行列Gとして出力 (引数渡し版)
	dispmat(G);
	dispmat(shiftright(F));		// 行列Fを1列分右にシフトして出力 (戻り値渡し版)
	dispmat(shiftright(F, 2));	// 行列Fを2列分右にシフトして出力 (戻り値渡し版)
	constexpr auto Gx4 = shiftright(Fx, 4);				// コンパイル時に右にシフト
	dispmat(Gx4);
	
	// 連結関連の関数
	printf("\n★★★★★★★ 連結関連の関数\n");
	constexpr ArcsMat<2,3> Hx1 = {
		10, 11, 22,
		33, 44, 55
	};
	constexpr ArcsMat<4,3> Hx2 = {
		66, 77, 88,
		99, 10, 12,
		23, 34, 56,
		67, 78, 89
	};
	ArcsMat<6,3> H12;
	concatv(Hx1, Hx2, H12);		// 行列を縦に連結する (引数渡し版)
	dispmat(Hx1);
	dispmat(Hx2);
	dispmat(H12);
	dispmat(concatv(Hx2, Hx1));	// 行列を縦に連結する (戻り値渡し版)
	constexpr auto Hx12 = concatv(Hx1, Hx2);			// コンパイル時に連結
	dispmat(Hx12);
	constexpr ArcsMat<4,1> Hx3 = {
		1,
		2,
		3,
		4	
	};
	ArcsMat<4,4> H44;
	concath(Hx2, Hx3, H44);		// 行列を横に連結する (引数渡し版)
	dispmat(H44);
	dispmat(concath(Hx2, Hx3));	// 行列を横に連結する (戻り値渡し版)
	constexpr auto Hx44 = concath(Hx2, Hx3);			// コンパイル時に連結
	dispmat(Hx44);
	constexpr ArcsMat<2,1> Hx4 = {
		5,
		6
	};
	ArcsMat<6,4> H64;
	concat4(Hx1, Hx3, Hx2, Hx4, H64);		// 4つの行列を連結する (引数渡し版)
	dispmat(H64);
	dispmat(concat4(Hx1, Hx4, Hx2, Hx3));	// 4つの行列を連結する (戻り値渡し版)
	constexpr auto Hx64 = concat4(Hx1, Hx3, Hx2, Hx4);	// コンパイル時に連結
	dispmat(Hx64);
	
	// 対角要素操作関連の関数
	printf("\n★★★★★★★ 対角要素操作関連の関数\n");
	constexpr ArcsMat<3,1> jx1 = {1, 3, 5};
	dispmat(jx1);
	ArcsMat<3,3> Jy1;
	diag(jx1, Jy1);							// 対角要素に縦ベクトルの各要素を持つ行列を生成する (引数渡し版) 
	dispmat(Jy1);
	dispmat(diag(jx1));						// 対角要素に縦ベクトルの各要素を持つ行列を生成する (戻り値渡し版) 
	constexpr auto Jx1 = diag(jx1);						// コンパイル時に生成
	dispmat(Jx1);
	constexpr ArcsMat<3,3> Jx2 = {
		3, 1, 1,
		1, 9, 1,
		1, 1, 2
	};
	dispmat(Jx2);
	ArcsMat<3,1> jy2;
	getdiag(Jx2, jy2);						// 対角要素を縦ベクトルとして取得する (引数渡し版)
	dispmat(jy2);
	dispmat(getdiag(Jx2));					// 対角要素を縦ベクトルとして取得する (戻り値渡し版)
	constexpr auto jx2 = getdiag(Jx2);					// コンパイル時に取得
	dispmat(jx2);
	constexpr ArcsMat<3,2> Jx3 = {
		3, 5,
		5, 9,
		5, 5
	};
	dispmat(Jx3);
	dispmat(getdiag(Jx3));					// 縦長の場合の対角要素を無理やり取得
	dispmat(getdiag(tp(Jx3)));				// 横長の場合の対角要素を無理やり取得
	printf("trace(Jx2) = %f\n", trace(Jx2));// トレースを計算する (戻り値渡し版のみ)
	constexpr double trJx2 = trace(Jx2);				// コンパイル時にトレースを計算
	printf("trJx2 = %f\n", trJx2);
	printf("multdiag(Jx2) = %f\n", multdiag(Jx2));	// 対角要素の総積を計算する (戻り値渡し版のみ)
	constexpr double mdJx2 = multdiag(Jx2);				// コンパイル時に対角の総積を計算
	printf("mdJx2 = %f\n", mdJx2);
	
	// 最大・最小値探索関連の関数
	printf("\n★★★★★★★ 最大・最小値探索関連の関数\n");
	constexpr ArcsMat<7,1> k1 = {
		 2,
		 5,
		-9,
		 4,
		 7,
		-1,
		 6
	};
	dispmat(k1);
	dispmat(Fx);
	printf("max(k1) = %f\n", max(k1));			// 縦ベクトル要素の最大値を返す関数
	printf("max(tp(k1)) = %f\n", max(tp(k1)));	// 横ベクトル要素の最大値を返す関数
	printf("max(Fx) = %f\n", max(Fx));			// 行列要素の最大値を返す関数
	constexpr double kx1 = max(k1);						// コンパイル時に最大値を計算
	printf("kx1 = %f\n", kx1);
	size_t kx, ky;
	std::tie(kx, ky) = maxidx(Fx);				// 行列要素の最大値の要素番号を返す関数
	printf("maxidx(Fx) = %ld, %ld\n", kx, ky);
	constexpr auto kxy = maxidx(k1);					// コンパイル時に要素番号を計算
	printf("kxy = %ld, %ld\n", std::get<0>(kxy), std::get<1>(kxy));
	printf("min(k1) = %f\n", min(k1));			// 縦ベクトル要素の最小値を返す関数
	printf("min(tp(k1)) = %f\n", min(tp(k1)));	// 横ベクトル要素の最小値を返す関数
	printf("min(Fx) = %f\n", min(Fx));			// 行列要素の最小値を返す関数
	constexpr double kx2 = min(k1);						// コンパイル時に最小値を計算
	printf("kx2 = %f\n", kx2);
	std::tie(kx, ky) = minidx(Fx);				// 行列要素の最小値の要素番号を返す関数
	printf("minidx(Fx) = %ld, %ld\n", kx, ky);
	constexpr auto kxy2 = minidx(k1);					// コンパイル時に要素番号を計算
	printf("kxy2 = %ld, %ld\n", std::get<0>(kxy2), std::get<1>(kxy2));
	
	// 要素ごとの数学関数
	printf("\n★★★★★★★ 要素ごとの数学関数\n");
	constexpr ArcsMat<3,4> Ax1 = {
		       1.0,        2.0,        3.0,      4.0,
		    M_PI_4,     M_PI_2, 3.0*M_PI_4,	    M_PI,
		5.0*M_PI_4, 3.0*M_PI_2, 7.0*M_PI_4, 2.0*M_PI 
	};
	dispmatfmt(Ax1, "%8.3f");
	printf("sum(Ax1) = %f\n", sum(Ax1));	// 行列要素の総和 (戻り値渡し版のみ)
	ArcsMat<3,4> Y1;
	exp(Ax1, Y1);							// 行列要素の指数関数 (引数渡し版)
	dispmatfmt(exp(Ax1), "%8.3f");			// 行列要素の指数関数 (戻り値渡し版)
	constexpr auto Yx1 = exp(Ax1);						// コンパイル時に行列要素の指数関数を計算
	dispmatfmt(Yx1, "%8.3f");
	ArcsMat<3,4> Y2;
	log(Y1, Y2);							// 行列要素の対数関数(底e版) (引数渡し版)
	dispmatfmt(log(Y1), "%8.3f");			// 行列要素の対数関数(底e版) (戻り値渡し版)
	constexpr auto Yx2 = log(Yx1);						// コンパイル時に行列要素の対数関数(底e版)を計算
	dispmatfmt(Yx2, "%8.3f");
	log10(Ax1, Y1);							// 行列要素の対数関数(底10版) (引数渡し版)
	dispmatfmt(log10(Ax1), "%8.3f");		// 行列要素の対数関数(底10版) (戻り値渡し版)
	constexpr auto Yx3 = log10(Ax1);					// コンパイル時に行列要素の対数関数(底10版)を計算
	dispmatfmt(Yx3, "%8.3f");
	sin(Ax1, Y1);							// 行列要素の正弦関数 (引数渡し版)
	dispmatfmt(sin(Ax1), "%8.3f");			// 行列要素の正弦関数 (戻り値渡し版)
	constexpr auto Yx4 = sin(Ax1);						// コンパイル時に行列要素の正弦関数を計算
	dispmatfmt(Yx4, "%8.3f");
	cos(Ax1, Y1);							// 行列要素の余弦関数 (引数渡し版)
	dispmatfmt(cos(Ax1), "%8.3f");			// 行列要素の余弦関数 (戻り値渡し版)
	constexpr auto Yx5 = cos(Ax1);						// コンパイル時に行列要素の余弦関数を計算
	dispmatfmt(Yx5, "%8.3f");
	tan(Ax1, Y1);							// 行列要素の正接関数 (引数渡し版)
	dispmatfmt(tan(Ax1), "% 12.3e");		// 行列要素の正接関数 (戻り値渡し版)
	constexpr auto Yx6 = tan(Ax1);						// コンパイル時に行列要素の正接関数を計算
	dispmatfmt(Yx6, "% 12.3e");
	sqrt(Ax1, Y1);							// 行列要素の平方根 (引数渡し版)
	dispmatfmt(sqrt(Ax1), "%8.3f");			// 行列要素の平方根 (戻り値渡し版)
	constexpr auto Yx7 = sqrt(Ax1);						// コンパイル時に行列要素の平方根を計算
	dispmatfmt(Yx7, "%8.3f");
	abs(Yx5, Y1);							// 行列要素の絶対値 (引数渡し版)
	dispmatfmt(abs(Yx5), "% 12.3f");		// 行列要素の絶対値 (戻り値渡し版)
	constexpr auto Yx8 = abs(Yx5);						// コンパイル時に行列要素の絶対値を計算
	dispmatfmt(Yx8, "% 12.3f");

	// 複素数関連の関数
	printf("\n★★★★★★★ 複素数関連の関数\n");
	constexpr ArcsMat<3,5,std::complex<double>> Acmpx1 = {
		-2, -1,  0,  1,  2,
		-4, -3,  0,  3,  4,
		-1, -2, -1, -2, -1
	};
	dispmat(Acmpx1);
	ArcsMat<3,5,std::complex<double>> Y9;
	sqrt(Acmpx1, Y9);						// 負の行列要素の平方根 (引数渡し版)
	dispmatfmt(sqrt(Acmpx1), "%6.3f");		// 負の行列要素の平方根 (戻り値渡し版)
	//constexpr auto Yx9 = sqrt(Acmpx1);				// コンパイル時に行列要素の平方根を計算 驚異の非対応エラー
	//dispmatfmt(Yx9, "%8.3f");
	constexpr ArcsMat<2,3,std::complex<double>> Acmpx2 = {
		std::complex( 1.0, 1.0), std::complex( 3.0, 4.0), std::complex( 3.0,-4.0),
		std::complex(-1.0, 1.0), std::complex(-3.0, 4.0), std::complex(-3.0,-4.0)
	};
	dispmat(Acmpx2);
	dispmat(abs(Acmpx2));

	// ノルム関連の関数
	printf("\n★★★★★★★ ノルム関連の関数\n");
	printf("norm<∞>(Ax1) = %f\n", norm<NormType::AMT_INFINITY>(Ax1));	// 無限大ノルムを計算する (戻り値渡し版のみ)
	constexpr double normAx1 = norm<NormType::AMT_INFINITY>(Ax1);				// コンパイル時に無限大ノルムを計算
	printf("norm<∞>(Ax1) = %f\n", normAx1);
/*
	// 三角行列操作関連の関数
	printf("\n★★★★★★★ 三角行列操作関連の関数\n");
	ArcsMat<7,6> Fx7;
	gettriup(Fx, Fx7);						// 上三角行列を切り出す (引数渡し版)
	dispmatfmt(Fx7, "%3.0f");
	dispmatfmt(gettriup(Fx, 3), "%3.0f");	// 上三角行列を切り出す (戻り値渡し版)
	constexpr auto Fxtri = gettriup(Fx, 2);				// コンパイル時に上三角行列を切り出す 
	dispmatfmt(Fxtri, "%3.0f");
	Fx7.FillAllZero();
	gettrilo(Fx, Fx7);						// 下三角行列を切り出す (引数渡し版)
	dispmatfmt(Fx7, "%3.0f");
	dispmatfmt(gettrilo(Fx, 3), "%3.0f");	// 下三角行列を切り出す (戻り値渡し版)
	constexpr auto Fxtri2 = gettrilo(Fx, 2);			// コンパイル時に下三角行列を切り出す 
	dispmatfmt(Fxtri2, "%3.0f");

	// LU分解関連の関数
	printf("\n★★★★★★★ LU分解関連の関数\n");
	ArcsMat<3,3> L, U, P;
	A.Set(
		10, -7,  0,
		-3,  2,  6,
		 5, -1,  5
	);
	dispmatfmt(A, "%3.0f");
	LUP(A, L, U, P);				// LU分解の結果と置換行列を計算 (引数渡し版)
	dispmatfmt(L, "%6.2f");
	dispmatfmt(U, "%6.2f");
	dispmatfmt(P, "%3.0f");
	dispmatfmt(tp(P)*L*U, "%3.0f");	// もとに戻るかチェック
	std::tie(L, U, P) = LUP(A);		// LU分解の結果と置換行列を計算 (タプル返し版)
	dispmatfmt(L, "%6.2f");
	dispmatfmt(U, "%6.2f");
	dispmatfmt(P, "%3.0f");
	LU(A, L, U);					// LU分解の結果のみを計算 (引数渡し版)
	dispmatfmt(L, "%6.2f");
	dispmatfmt(U, "%6.2f");
	dispmatfmt(L*U, "%3.0f");		// もとに戻るかチェック
	std::tie(L, U) = LU(A);			// LU分解の結果のみを計算 (タプル返し版)
	dispmatfmt(L, "%6.2f");
	dispmatfmt(U, "%6.2f");
	dispmatfmt(L*U, "%3.0f");		// もとに戻るかチェック
	constexpr auto LxUx = LU(Ax);						// コンパイル時にLU分解を計算
	constexpr auto Lx = std::get<0>(LxUx);				// コンパイル時に計算した下三角を抽出
	constexpr auto Ux = std::get<1>(LxUx);				// コンパイル時に計算した上三角を抽出
	dispmatfmt(Ax, "%3.0f");
	dispmatfmt(Lx, "%6.2f");
	dispmatfmt(Ux, "%6.2f");
	dispmatfmt(Lx*Ux, "%3.0f");		// もとに戻るかチェック
*/
	/*
	// 行列演算補助関連の関数のテスト
	printf("\n★★★★★★★ 行列演算補助関連の関数のテスト\n");
	printf("nonzeroele(I) = %ld\n", nonzeroele(I));
	printf("rank(I) = %ld\n", rank(I));
	
	// ノルム演算関連のテスト
	printf("\n★★★★★★★ ノルム演算関連のテスト\n");
	printf("infnorm(A) = %f\n", infnorm(A));
	printf("euclidnorm(v) = %f\n", euclidnorm(v));
	
	// Cholesky分解(LDL^T版)のテスト
	printf("\n★★★★★★★ Cholesky分解(LDL^T版)のテスト\n");
	Matrix<3,3> Ach = {
		  4,  12, -16,
		 12,  37, -43,
		-16, -43,  98
	};
	Matrix<3,3> Lch, Dch;
	Cholesky(Ach, Lch, Dch);
	PrintMat(Ach);
	PrintMat(Lch);
	PrintMat(Dch);
	PrintMat(Lch*Dch*tp(Lch));
	
	// Cholesky分解(LL^T版)のテスト
	printf("\n★★★★★★★ Cholesky分解(LL^T版)のテスト\n");
	Cholesky(Ach, Lch);
	PrintMat(Ach);
	PrintMat(Lch);
	PrintMat(Lch*tp(Lch));
	
	// QR分解のテスト1
	printf("\n★★★★★★★ QR分解(実数版)のテスト1\n");
	Matrix<3,3> Aqr = {
		2, -2, 18,
		2,  1,  0,
		1,  2,  0
	};
	Matrix<3,3> Qqr, Rqr;
	QR(Aqr, Qqr, Rqr);
	PrintMat(Aqr);
	PrintMatrix(Qqr, "% 8.3f");
	PrintMatrix(Rqr, "% 8.3f");
	PrintMatrix(Qqr*tp(Qqr), "% 7.3f");	// Qが直交行列かチェック
	PrintMat(Qqr*Rqr);					// 元に戻るかチェック
	
	// QR分解のテスト2
	printf("\n★★★★★★★ QR分解(実数版)のテスト2\n");
	Aqr.Set(
		12, -51,   4,
		 6, 167, -68,
		-4,  24, -41
	);
	QR(Aqr, Qqr, Rqr);
	PrintMat(Aqr);
	PrintMatrix(Qqr, "% 8.3f");
	PrintMatrix(Rqr, "% 8.3f");
	PrintMatrix(Qqr*tp(Qqr), "% 7.3f");	// Qが直交行列かチェック
	PrintMat(Qqr*Rqr);					// 元に戻るかチェック
	
	// QR分解のテスト3
	printf("\n★★★★★★★ QR分解(実数版)のテスト3\n");
	Matrix<4,3> Aqr3 = {
		12, -51,   4, 39,
		 6, 167, -68, 22,
		-4,  24, -41, 11
	};
	Matrix<3,3> Qqr3;
	Matrix<4,3> Rqr3;
	QR(Aqr3, Qqr3, Rqr3);
	PrintMat(Aqr3);
	PrintMatrix(Qqr3, "% 8.3f");
	PrintMatrix(Rqr3, "% 8.3f");
	PrintMatrix(Qqr3*tp(Qqr3), "% 7.3f");	// Qが直交行列かチェック
	PrintMat(Qqr3*Rqr3);					// 元に戻るかチェック
	
	// QR分解のテスト4
	printf("\n★★★★★★★ QR分解(実数版)のテスト4\n");
	Matrix<3,4> Aqr4 = tp(Aqr3);
	Matrix<4,4> Qqr4;
	Matrix<3,4> Rqr4;
	QR(Aqr4, Qqr4, Rqr4);
	PrintMat(Aqr4);
	PrintMatrix(Qqr4, "% 8.3f");
	PrintMatrix(Rqr4, "% 8.3f");
	PrintMatrix(Qqr4*tp(Qqr4), "% 7.3f");	// Qが直交行列かチェック
	PrintMat(Qqr4*Rqr4);					// 元に戻るかチェック
	
	// SVD特異値分解のテスト1(縦長行列の場合)
	printf("\n★★★★★★★ SVD特異値分解のテスト1(縦長行列の場合)\n");
	Matrix<2,4> As = {
		1, 2,
		3, 4,
		5, 6,
		7, 8
	};
	Matrix<4,4> Us;
	Matrix<2,4> Ss;
	Matrix<2,2> Vs;
	SVD(As, Us, Ss, Vs);
	PrintMat(As);
	PrintMat(Us);
	PrintMat(Ss);
	PrintMat(Vs);
	PrintMat(Us*Ss*tp(Vs));	// 元に戻るかチェック
	
	// SVD特異値分解のテスト2(横長行列の場合)
	printf("\n★★★★★★★ SVD特異値分解のテスト2(横長行列の場合)\n");
	Matrix<4,2> As2 = tp(As);
	Matrix<2,2> Us2;
	Matrix<4,2> Ss2;
	Matrix<4,4> Vs2;
	SVD(As2, Us2, Ss2, Vs2);
	PrintMat(As2);
	PrintMat(Us2);
	PrintMat(Ss2);
	PrintMat(Vs2);
	PrintMat(Us2*Ss2*tp(Vs2));	// 元に戻るかチェック
	
	// SVD特異値分解のテスト3(ランク落ちの場合)
	printf("\n★★★★★★★ SVD特異値分解のテスト3(ランク落ちの場合)\n");
	Matrix<3,3> As3 = {
		 2,  0,  2,
		 0,  1,  0,
		 0,  0,  0
	};
	Matrix<3,3> Us3;
	Matrix<3,3> Ss3;
	Matrix<3,3> Vs3;
	SVD(As3, Us3, Ss3, Vs3);
	PrintMat(As3);
	PrintMat(Us3);
	PrintMat(Ss3);
	PrintMat(Vs3);
	PrintMat(Us3*Ss3*tp(Vs3));	// 元に戻るかチェック
	printf("rank(As3) = %ld\n", rank(As3));
	
	// SVD特異値分解のテスト4(符号修正が必要な場合)
	printf("\n★★★★★★★ SVD特異値分解のテスト4(符号修正が必要な場合)\n");
	Matrix<3,3> As4 = {
		 1,  1,  3,
		-5,  6, -3, 
		 7, -2,  9
	};
	Matrix<3,3> Us4;
	Matrix<3,3> Ss4;
	Matrix<3,3> Vs4;
	SVD(As4, Us4, Ss4, Vs4);
	PrintMat(As4);
	PrintMat(Us4);
	PrintMat(Ss4);
	PrintMat(Vs4);
	PrintMat(Us4*Ss4*tp(Vs4));	// 元に戻るかチェック
	
	// Schur分解のテスト1(実数固有値の場合)
	printf("\n★★★★★★★ Schur分解のテスト1(実数固有値の場合)\n");
	Matrix<3,3> Asr1 = {
		 4, -8,  1,
		 1, -5,  1,
		-9,  8, -6
	};
	Matrix<3,3> Qsr, Usr;
	Schur(Asr1, Qsr, Usr);
	PrintMat(Asr1);
	PrintMatrix(Qsr, "%6.3f");
	PrintMatrix(Usr, "%6.3f");
	PrintMat(Qsr*Usr*inv(Qsr));
	
	// Schur分解のテスト2(複素数固有値の場合)
	printf("\n★★★★★★★ Schur分解のテスト2(複素数固有値の場合)\n");
	Matrix<3,3> Asr2 = {
		 1,  2,  3,
		-5,  9, -1,
		 2,  6,  8
	};
	Schur(Asr2, Qsr, Usr);
	PrintMat(Asr2);
	PrintMatrix(Qsr, "%6.3f");
	PrintMatrix(Usr, "%6.3f");
	PrintMat(Qsr*Usr*inv(Qsr));
	
	// 連立方程式の球解テスト
	printf("\n★★★★★★★ 連立方程式の球解テスト\n");
	Matrix<1,3> b;
	b.Set(
		9,
		5,
		7
	);
	PrintMat(A);
	PrintMat(b);
	Matrix<1,3> xslv;
	solve(A, b, xslv);
	PrintMat(xslv);
	PrintMat(solve(A, b));
	
	// 上三角行列の連立方程式の球解テスト
	printf("\n★★★★★★★ 上三角行列の連立方程式の球解テスト\n");
	Matrix<3,3> Auptri = {
		1, 3, 6,
		0, 2, 7,
		0, 0,-4
	};
	solve_upper_tri(Auptri, b, xslv);
	PrintMat(Auptri);
	PrintMat(xslv);
	
	// 行列式と逆行列のテスト
	printf("\n★★★★★★★ 行列式と逆行列のテスト\n");
	printf("det(A) = %f\n", det(A));
	PrintMatrix(inv(A), "% 16.14e");
	PrintMatrix(inv_with_check(A), "% 16.14e");
	
	// 左上小行列の逆行列のテスト
	printf("\n★★★★★★★ 左上小行列の逆行列のテスト\n");
	Matrix<5,5> Ai5 = {
		1,  1,  1,  0,  0,
		2,  3, -2,  0,  0,
		3, -1,  1,  0,  0,
		0,  0,  0,  0,  0,
		0,  0,  0,  0,  0
	};
	PrintMat(Ai5);
	PrintMatrix(inv(Ai5, 3), "% 16.14e");
	
	// 上三角行列の逆行列のテスト
	printf("\n★★★★★★★ 上三角行列の逆行列のテスト\n");
	Matrix<3,3> Auptri_inv;
	inv_upper_tri(Auptri, Auptri_inv);
	PrintMat(Auptri_inv);
	
	// 上三角行列で左上小行列の逆行列のテスト
	printf("\n★★★★★★★ 上三角行列で左上小行列の逆行列のテスト\n");
	Matrix<4,4> Auptri2 = {
		1, 3, 6, 0,
		0, 2, 7, 0,
		0, 0,-4, 0,
		0, 0, 0, 0
	};
	Matrix<4,4> Auptri_inv2;
	PrintMat(Auptri2);
	inv_upper_tri(Auptri2, 3, Auptri_inv2);
	PrintMat(Auptri_inv2);
	
	// 疑似逆行列のテスト
	printf("\n★★★★★★★ 疑似逆行列のテスト\n");
	PrintMatrix(lpinv(D), "% 16.14e");
	PrintMatrix(rpinv(tp(D)), "% 16.14e");
	Matrix<1,2> Dpinv = {
		1,
		2
	};
	PrintMatrix(lpinv(Dpinv), "% 16.14e");
	PrintMatrix(rpinv(tp(Dpinv)), "% 16.14e");
	
	// 左上小行列の疑似逆行列のテスト
	printf("\n★★★★★★★ 左上小行列の疑似逆行列のテスト\n");
	Matrix<4,5> Dpinv45 = {
		1,  0,  0,  0,
		2,  0,  0,  0,
		0,  0,  0,  0,
		0,  0,  0,  0,
		0,  0,  0,  0,
	};
	PrintMat(Dpinv45);
	PrintMatrix(lpinv(Dpinv45, 1), "% 16.14e");
	PrintMatrix(rpinv(tp(Dpinv45), 1), "% 16.14e");
	Dpinv45.Set(
		1,  1,  0,  0,
		2,  3,  0,  0,
		3, -1,  0,  0,
		0,  0,  0,  0,
		0,  0,  0,  0
	);
	PrintMat(Dpinv45);
	PrintMatrix(lpinv(Dpinv45, 2), "% 16.14e");
	PrintMatrix(rpinv(tp(Dpinv45), 2), "% 16.14e");
	
	// 行列指数関数のテスト
	printf("\n★★★★★★★ 行列指数関数のテスト\n");
	Matrix<3,3> Y = expm(A, 6);
	PrintMat(A);
	PrintMatrix(Y,"% 16.14e");
	printf("det(Y)     = % 16.14e\n", det(Y));
	printf("exp(tr(A)) = % 16.14e\n\n", exp(tr(A)));	// 公式通りに一致！
	PrintMatrix(integral_expm(A,100e-6,10,6),"% 16.14e");
	
	// 特に意味のないリッカチ方程式のテスト
	printf("\n★★★★★★★ 特に意味のないリッカチ方程式のテスト\n");
	Matrix<3,3> P;
	P.Set(
		1,  2, -5,
		2,  1, -1,
		7, -6,  9
	);
	PrintMat(P);
	auto Q = P*2.71828;
	auto DP = P*3.14159;
	auto Z = P*A + tp(A)*P - P*B*tp(B)*P + Q + DP;
	PrintMatrix(Z,"% 10.3f");
	
	// float型のテスト
	printf("★★★ float型のテスト\n");
	Matrix<3,3,float> Af = {
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	};
	PrintMat(Af);
	auto Yf = expm(Af, 6);
	PrintMatrix(Yf, "%16.14e");
	printf("det(Y)     = %16.14e\n", det(Yf));
	printf("exp(tr(A)) = %16.14e\n\n", exp(tr(Af)));	// 公式通りに一致！
	PrintMatrix(integral_expm(Af,100e-6,10,6),"% 16.14e");
	
	// int型のテスト
	printf("★★★ int型のテスト\n");
	Matrix<2,2,int> Ai = {
		1, 2,
		3, 4
	};
	Matrix<2,2,int> Bi = {
		5, 6,
		7, 8
	};
	PrintMat(Ai);
	PrintMat(Bi);
	PrintMat(Ai*Bi);
	
	// long型のテスト
	printf("★★★ long型のテスト\n");
	Matrix<2,2,long> Al = {
		1, 2,
		3, 4
	};
	Matrix<2,2,long> Bl = {
		5, 6,
		7, 8
	};
	PrintMat(Al);
	PrintMat(Bl);
	PrintMat(Al*Bl);
	
	// 複素数型のテスト
	printf("★★★ 複素数型のテスト\n");
	Matrix<2,2,std::complex<double>> Ac = {
		std::complex(1.0,2.0), std::complex(2.0,3.0),
		std::complex(3.0,4.0), std::complex(4.0,5.0)
	};
	Matrix<2,2,std::complex<double>> Bc = {
		std::complex(1.0,2.0), std::complex(2.0,3.0),
		std::complex(3.0,4.0), std::complex(4.0,5.0)
	};
	PrintMat(Ac);
	PrintMat(Bc);
	PrintMat(Ac*Bc);
	PrintMat(Ac*Bc - std::complex(0.0,29.0));
	PrintMat(reale(Ac));	// 実数部
	PrintMat(image(Ac));	// 虚数部
	PrintMat(mage(Ac));		// 大きさ
	PrintMat(arge(Ac));		// 偏角
	PrintMat(conje(Ac));	// 複素共役
	PrintMat((Matrix<3,3,std::complex<double>>::eye()));	// 単位行列
	
	// 負の平方根のテスト
	printf("★★★ 負の平方根のテスト\n");
	Matrix<3,3,std::complex<double>> Acomp = {
		std::complex(1.0,0.0), std::complex( 1.0,0.0), std::complex( 1.0,0.0),
		std::complex(2.0,0.0), std::complex( 3.0,0.0), std::complex(-2.0,0.0),
		std::complex(3.0,0.0), std::complex(-1.0,0.0), std::complex( 1.0,0.0),
	};
	PrintMat(Acomp);
	PrintMatrix(sqrte(A), "% 6.3f");	// double型だとnanになってしまうが，
	PrintMatrix(sqrte(Acomp), "% 6.3f");// std::complexだとちゃんと計算できる
	
	// エルミート転置のテスト
	printf("★★★ エルミート転置のテスト\n");
	Matrix<2,3,std::complex<double>> Acomp2 = {
		std::complex(1.0, 2.0), std::complex(  3.0, -4.0),
		std::complex(5.0,-6.0), std::complex(  7.0,  8.0),
		std::complex(9.0,10.0), std::complex(-11.0,-12.0)
	};
	PrintMat(Acomp2);
	PrintMat(Htp(Acomp2));
	
	// 複素数LU分解のテスト
	printf("★★★ 複素数LU分解のテスト\n");
	Matrix<3,3,std::complex<double>> Acomp3 = {
		std::complex( 4.0, 6.0), std::complex( 1.0,-3.0), std::complex( 5.0, 2.0),
		std::complex( 8.0,-5.0), std::complex(-7.0,-6.0), std::complex( 7.0,-1.0),
		std::complex( 9.0, 9.0), std::complex(-7.0,-5.0), std::complex(-5.0,-3.0)
	};
	Matrix<3,3,std::complex<double>> Lcomp, Ucomp;
	Matrix<1,3,int> vcomp;
	LU(Acomp3, Lcomp, Ucomp, vcomp);
	PrintMat(Acomp3);
	PrintMat(Lcomp);
	PrintMat(Ucomp);
	PrintMat(reorderrow(Lcomp*Ucomp, vcomp));
	
	// 複素数逆行列のテスト
	printf("★★★ 複素数逆行列のテスト\n");
	PrintMat(inv(Acomp3));
	PrintMatrix(inv(Acomp3)*Acomp3, "% 7.3f");
	
	// 複素数QR分解のテスト
	printf("★★★ 複素数QR分解のテスト\n");
	Matrix<3,3,std::complex<double>> Qcqr, Rcqr;
	QR(Acomp3, Qcqr, Rcqr);
	PrintMatrix(Qcqr, "% 7.3f");
	PrintMatrix(Rcqr, "% 7.3f");
	PrintMatrix(Qcqr*Htp(Qcqr), "% 7.3f");	// Qが直交行列かチェック
	PrintMat(Qcqr*Rcqr);					// 元に戻るかチェック
	
	// 固有値計算のテスト1(実数固有値の場合)
	printf("★★★ 固有値計算のテスト1(実数固有値の場合)\n");
	Matrix<3,3> Aeig = {
		-3, -4,  2,
		-7,  1, -5,
		 6, -7,  3
	};
	PrintMat(Aeig);
	PrintMat(eigen(Aeig));
	PrintMat(eigenvec(Aeig));
	
	// 固有値計算のテスト2(複素数固有値の場合)
	printf("★★★ 固有値計算のテスト2(複素数固有値の場合)\n");
	Aeig.Set(
		10,  -8,   5,
		-8,   9,   6,
		-1, -10,   7
	);
	PrintMat(Aeig);
	PrintMat(eigen(Aeig));
	PrintMat(eigenvec(Aeig));
	
	// 2乗のコンパイル時定数演算のテスト
	printf("\n★★★★★★★ 2乗のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cx = {
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	};
	constexpr Matrix<3,3> Cxsq = Cx^2;			// 2乗のコンパイル時計算
	PrintMat(Cxsq);
	
	// 逆行列のコンパイル時定数演算のテスト
	printf("\n★★★★★★★ 逆行列のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cxinv = inv(Cx);		// 逆行列のコンパイル時計算
	PrintMat(Cxinv);
	
	// 状態遷移行列のコンパイル時定数演算のテスト
	printf("\n★★★★★★★ 状態遷移行列のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cxexp = expm(Cx, 6);	// 状態遷移行列のコンパイル時計算
	PrintMatrix(Cxexp, "% 16.14e");
	printf("det(Y)     = % 16.14e\n", det(Cxexp));
	printf("exp(tr(A)) = % 16.14e\n\n", exp(tr(Cx)));	// 公式通りに一致！
	
	// 状態遷移行列の定積分のコンパイル時定数演算のテスト
	printf("\n★★★★★★★ 状態遷移行列の定積分のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cxexpint = integral_expm(Cx, 100e-6, 10, 6);	// 状態遷移行列の定積分のコンパイル時計算
	PrintMatrix(Cxexpint, "% 16.14e");
	
	// SVD特異値分解のコンパイル時定数演算のテスト
	printf("\n★★★★★★★ SVD特異値分解のコンパイル時定数演算のテスト\n");
	constexpr Matrix<2,4> Axsvd = {
		1, 2,
		3, 4,
		5, 6,
		7, 8
	};
	constexpr auto USVx = SVD(Axsvd);
	constexpr Matrix<4,4> Uxsvd = std::get<0>(USVx);
	constexpr Matrix<2,4> Sxsvd = std::get<1>(USVx);
	constexpr Matrix<2,2> Vxsvd = std::get<2>(USVx);
	PrintMat(Uxsvd);
	PrintMat(Sxsvd);
	PrintMat(Vxsvd);
	PrintMat(Uxsvd*Sxsvd*tp(Vxsvd));	// 元に戻るかチェック
	
	// 行列のランクのコンパイル時定数演算のテスト
	printf("\n★★★★★★★ 行列のランクのコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Ark = {
		 2,  0,  2,
		 0,  1,  0,
		 0,  0,  0
	};
	constexpr size_t RankOfArk = rank(Ark);
	PrintMat(Ark);
	printf("rank(Ark) = %zu\n", RankOfArk);
	static_assert(rank(Ark) == 2);
	
	// クロネッカー積のテスト
	printf("\n★★★★★★★ クロネッカー積のテスト\n");
	constexpr Matrix<2,2> Ukl = {
		1, 2,
		3, 4
	};
	constexpr Matrix<2,2> Ukr = {
		0, 5,
		6, 7
	};
	PrintMat(Kronecker(Ukl, Ukr));
	PrintMat(A);
	PrintMat(Axsvd);
	PrintMat(Kronecker(A, Axsvd));
	PrintMat(Kronecker(Axsvd, A));
	
	// vec作用素のテスト
	printf("\n★★★★★★★ vec作用素のテスト\n");
	PrintMat(vec(A));
	PrintMat(vec(Axsvd));
	Matrix<1,9>::vecinv(vec(A), A);			// 引数で返す版
	A = Matrix<1,9>::vecinv<3,3>(vec(A));	// 戻り値で返す版
	PrintMat(A);
	auto Avec = Matrix<1,8>::vecinv<2,4>(vec(Axsvd));	// 戻り値で返す版
	PrintMat(Avec);
	*/
	return EXIT_SUCCESS;	// 正常終了
}

