//! @file OfflineFunction.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2024/06/26
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
	ArcsMat<3,3> A = {			// 宣言と同時に値をセットする場合
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	};
	ArcsMat<3,3> B;				// 宣言したあとに、
	B.Set(						// 値をセットする場合
		5,  1, -3,
		3, -1,  7,
		4,  2, -5
	);
	auto C = A;					// 宣言と同時に既にある行列を代入する場合
	ArcsMat<2,3> Pi(M_PI);		// 宣言と同時に設定値で全て埋める場合
	ArcsMat<2,2,int> Aint = {	// int型の行列を定義する場合
		 1, 3,
		-5, 7
	};
	//ArcsMat<2,2,std::string> Astr;	// 数値型以外だとエラー！ (アンコメントするとAssertion Failed)

	// 行列の表示
	printf("\n★★★★★★★ 行列の表示\n");
	dispsize(A);			// 行列Aのサイズを表示
	disp(A);				// 行列Aを表示
	disp(B);				// 行列Bを表示
	dispf(C, "% 6.4f");		// 表示の書式指定をして表示する場合
	dispf(Pi, "% 16.15f");	// 行列Piを書式指定して表示
	disp(Aint);				// 整数行列Aintを表示
	A.DispAddress();		// 行列Aのメモリアドレスを表示

	// 行列のサイズと要素情報の取得
	printf("H = %zu, W = %zu\n", A.GetHeight(), A.GetWidth());	// 行列Aの高さと幅を表示
	printf("Non-Zero = %zu\n", A.GetNumOfNonZero());			// 行列Aの非ゼロ要素の数を表示

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
		disp(J[i]);
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
		disp(K.at(i));
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
	disp(Alpha);
	disp(Ax);
	
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
	disp(alpha);
	printf("alpha[3] = % g\n", alpha[3]);	// 縦ベクトルの3番目を表示
	alpha[3] = 6;			// 縦ベクトルの3番目に 6 を書き込む
	disp(alpha);			// 縦ベクトルの3番目が変わっている
	
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
	disp(B);	// 変わっている
	B(2,3) = 7;	// 2行3列目に 7 を書き込む
	disp(B);	// 元に戻っている
	
	// 行列の()オペレータによる要素アクセス(サイズチェック可能版)
	printf("\n★★★★★★★ 行列の()オペレータによる要素アクセス(サイズチェック可能版)\n");
	printf("B(2,3,true) = % g\n", B(2,3,true));		// 何も問題ない場合 (trueを付けるとサイズチェック有り)
	//printf("B(2,4,true) = % g\n", B(2,4,true));	// ←サイズエラー！
	//printf("B(2,4,true) = % g\n", B(2,4,false));	// ←運が悪いとSIGSEGV (falseにするとサイズチェック無し)
	B(2,3,true) = 0;	// 2行3列目に 0 を書き込む
	disp(B);			// 変わっている
	B(2,3,true) = 7;	// 2行3列目に 7 を書き込む
	disp(B);			// 元に戻っている
	//B(2,4,true) = 7;	// ←サイズエラー！
	//B(2,4,false) = 7;	// ←運が悪いとSIGSEGV
	
	// 代入演算子
	printf("\n★★★★★★★ 代入演算子\n");
	C = B;
	disp(C);
	//C = Alpha;	// ←サイズエラー！
	
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
	//C = A/B;			// ←この演算子の使い方は使用禁止！
	//C = 1/A;			// ←この演算子の使い方は使用禁止！
	
	// すべての要素を埋める
	printf("\n★★★★★★★ すべての要素を埋める\n");
	C.FillAll(3.9);	// すべての要素を3.9で埋める
	disp(C);
	C.FillAllZero();// すべての要素をゼロで埋める
	disp(C);
	C = B;
	
	// 1次元std::array配列との相互変換
	printf("\n★★★★★★★ 1次元std::array配列との相互変換\n");
	std::array<double, 10> v10array = {100,200,300,400,500,600,700,800,900,1000};
	ArcsMat<10,1> v10;
	v10.LoadArray(v10array);	// std::arrayから読み込む
	disp(v10);
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
	disp(D);
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
	disp(C);
	ArcsMat<2,1> v2;
	C.GetVerticalVec(v2,2,3);			// 2行3列目を先頭として高さ2の縦ベクトルを抜き出す (引数渡し版)
	disp(v2);
	disp(C.GetVerticalVec<2>(2,3));		// 2行3列目を先頭として高さ2の縦ベクトルを抜き出す (戻り値渡し版)
	ArcsMat<1,2> w2;
	C.GetHorizontalVec(w2,2,2);			// 2行2列目を先頭として幅2の横ベクトルを抜き出す (引数渡し版)
	disp(w2);	
	disp(C.GetHorizontalVec<2>(2,2));	// 2行2列目を先頭として幅2の横ベクトルを抜き出す (戻り値渡し版)
	//disp(C.GetVerticalVec<2>(3,3));	// ←はみ出るエラー！(アンコメントするとAssertion Failed)
	//disp(C.GetHorizontalVec<2>(2,3));	// ←はみ出るエラー！(アンコメントするとAssertion Failed)
	constexpr auto v2x = Ax.GetVerticalVec<2>(1,3);		// コンパイル時に抜き出す
	constexpr auto w2x = Ax.GetHorizontalVec<2>(3,1);	// コンパイル時に抜き出す
	disp(v2x);
	disp(w2x);

	// ゼロの取り扱い関連の関数
	printf("\n★★★★★★★ ゼロの取り扱い関連の関数\n");
	const ArcsMat<3,3> Azero1 = {
		 ArcsMat<3,3>::EPSILON, -ArcsMat<3,3>::EPSILON, ArcsMat<3,3>::EPSILON,
		-ArcsMat<3,3>::EPSILON,  ArcsMat<3,3>::EPSILON, ArcsMat<3,3>::EPSILON,
		 ArcsMat<3,3>::EPSILON, -ArcsMat<3,3>::EPSILON, ArcsMat<3,3>::EPSILON
	};
	dispf(Azero1, "% 16.15e");
	auto Azero2 = Azero1;
	Azero2.Zeroing();			// ゼロに近い要素を完全にゼロにする
	dispf(Azero2, "% 16.15e");
	Azero2.Zeroing(1e-12);		// ゼロに近い要素を完全にゼロにする (許容誤差指定版)
	dispf(Azero2, "% 16.15e");
	Azero2 = Azero1;
	Azero2.ZeroingTriLo(1e-12);	// 下三角(主対角除く)に限定して、ゼロに近い要素を完全にゼロにする (許容誤差指定版)
	dispf(Azero2, "% 16.15e");
	const ArcsMat<2,2,std::complex<double>> Azero3 = {
		std::complex(ArcsMat<3,3>::EPSILON, 1.0), std::complex(1.0, ArcsMat<3,3>::EPSILON),
		                  ArcsMat<3,3>::EPSLCOMP, ArcsMat<3,3>::EPSLCOMP
	};
	dispf(Azero3, "% 16.15e");
	auto Azero4 = Azero3;
	Azero4.Zeroing(1e-12);		// ゼロに近い要素を完全にゼロにする (複素数の場合)
	dispf(Azero4, "% 16.15e");
	Azero4 = Azero3;
	Azero4.ZeroingTriLo(1e-12);	// 下三角(主対角除く)に限定して、ゼロに近い要素を完全にゼロにする (複素数の場合)
	dispf(Azero4, "% 16.15e");

	// ベクトルの埋め込み
	printf("\n★★★★★★★ ベクトルの埋め込み\n");
	v2.Set(
		39,
		93
	);
	C.SetVerticalVec(v2,2,3);		// 2行3列目を先頭として縦ベクトルを埋め込む
	disp(C);
	w2.Set(39, 93);
	C.SetHorizontalVec(w2,2,1);		// 2行1列目を先頭として横ベクトルを埋め込む
	disp(C);
	//C.SetVerticalVec(v2,3,1);		// ←はみ出るエラー！(アンコメントするとAssertion Failed)
	//C.SetHorizontalVec(w2,1,3);	// ←はみ出るエラー！(アンコメントするとAssertion Failed)

	// 零行列と1行列と単位行列
	printf("\n★★★★★★★ 零行列と1行列と単位行列\n");
	auto O = zeros<5,5>();	// 零行列
	auto l = ones<5,5>();	// 1行列
	auto I = eye<5,5>();	// 単位行列
	disp(O);
	disp(l);
	disp(I);
	constexpr auto Ox = zeros<3,5>();	// コンパイル時に零行列を生成
	constexpr auto lx = ones<3,5>();	// コンパイル時に1行列を生成
	constexpr auto Ix = eye<3,3>();		// コンパイル時に単位行列を生成
	disp(Ox);
	disp(lx);
	disp(Ix);
	
	// 単調増加ベクトル
	printf("\n★★★★★★★ 単調増加ベクトル\n");
	auto r = ramp<5,1>();				// 単調増加ベクトル
	constexpr auto rx = ramp<5,1>();	// コンパイル時に単調増加ベクトルを生成
	disp(r);
	disp(rx);
	
	// 列操作関連の関数
	printf("\n★★★★★★★ 列操作関連の関数\n");
	ArcsMat<3,1> vd;
	getcolumn(D, vd, 2);			// 行列Dの2列目を縦ベクトルとして抽出 (引数渡し版)
	disp(vd);
	disp(getcolumn(D, 1));			// 行列Dの1列目を縦ベクトルとして抽出 (戻り値渡し版)
	constexpr auto vax = getcolumn(Ax, 2);				// コンパイル時に縦ベクトルを抽出
	disp(vax);
	setcolumn(D, vax, 1);			// 行列Dの1列目に縦ベクトルを埋め込む (引数渡し版)
	disp(D);
	disp(setcolumn(vd, 2, D));		// 行列Dの2列目に縦ベクトルを埋め込む (戻り値渡し版)
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
	disp(Ex);
	disp(vex);
	constexpr auto Ex2 = setcolumn(vex, 3, Ex);			// コンパイル時に縦ベクトルを埋め込む
	disp(Ex2);
	ArcsMat<3,4> E = Ex;
	disp(E);
	swapcolumn(E, 1, 3);				// 行列Eの1列目と3列目を入れ替える (引数渡し版)
	disp(E);
	disp(swapcolumn(1, 3, E));			// 行列Eの1列目と3列目を入れ替える (戻り値渡し版)
	constexpr auto Ex3 = swapcolumn(1, 3, Ex2);			// コンパイル時に列を入れ替える
	disp(Ex3);
	fillcolumn(E, 39, 4, 2, 3);			// 行列Eの4列目2行目～3行目を39で埋める (引数渡し版)
	disp(E);
	disp(fillcolumn(1, 4, 2, 3, E));	// 行列Eの4列目2行目～3行目を1で埋める (戻り値渡し版)
	constexpr auto Ex4 = fillcolumn(39, 3, 1, 2, Ex3);	// コンパイル時に指定列を値で埋める
	disp(Ex4);
	ArcsMat<1,4,size_t> wor = {2, 1, 4, 3};
	disp(wor);
	E = ordercolumn(E, wor);		// 行列Eをworが昇順になるように列を並び替える (戻り値渡し版のみ)
	disp(E);
	E = ordercolumn(E, wor);		// 2回やると元に戻る
	disp(E);
	ordercolumn_and_vec(E, wor);	// 行列Eをworが昇順になるようにEとworの両方の列を並び替える (引数渡し版のみ)
	disp(E);
	disp(wor);
	constexpr ArcsMat<1,4,size_t> worx = {2, 1, 4, 3};
	constexpr auto Ex5 = ordercolumn(Ex, worx);			// コンパイル時に並び替える
	disp(Ex5);
	ArcsMat<1,4> e5;	
	sumcolumn(Ex5, e5);				// 各列の総和を計算する (引数渡し版)
	disp(e5);
	disp(sumcolumn(Ex5));			// 各列の総和を計算する (戻り値渡し版)
	constexpr auto ex5 = sumcolumn(Ex5);				// コンパイル時に総和を計算する
	disp(ex5);
	
	// 行操作関連の関数
	printf("\n★★★★★★★ 行操作関連の関数\n");
	ArcsMat<1,2> wd;
	getrow(D, wd, 2);				// 行列Dの2行目を横ベクトルとして抽出 (引数渡し版)
	disp(wd);
	disp(getrow(D, 3));				// 行列Dの3行目を横ベクトルとして抽出 (戻り値渡し版)
	constexpr auto wax = getrow(Ax, 3);					// コンパイル時に横ベクトルを抽出
	disp(wax);
	setrow(D, wd, 1);				// 行列Dの1行目に横ベクトルを埋め込む (引数渡し版)
	disp(D);
	disp(setrow(wd, 2, D));			// 行列Dの2行目に横ベクトルを埋め込む (戻り値渡し版)
	constexpr ArcsMat<4,3> Etx = {
		0, 1, 2,
		3, 4, 5,
		6, 7, 8,
		9, 0, 1
	};
	constexpr ArcsMat<1,3> wetx = {39, 39, 24};
	disp(Etx);
	disp(wetx);
	constexpr auto Etx2 = setrow(wetx, 3, Etx);			// コンパイル時に横ベクトルを埋め込む
	disp(Etx2);
	ArcsMat<4,3> Et = Etx;
	disp(Et);
	swaprow(Et, 1, 3);					// 行列Etの1行目と3行目を入れ替える (引数渡し版)
	disp(Et);
	disp(swaprow(1, 3, Et));			// 行列Etの1行目と3行目を入れ替える (戻り値渡し版)
	constexpr auto Etx3 = swaprow(1, 3, Etx2);			// コンパイル時に行を入れ替える
	disp(Etx3);
	fillrow(Et, 39, 4, 2, 3);			// 行列Etの4行目2列目～3列目を39で埋める (引数渡し版)
	disp(Et);
	disp(fillrow(1, 4, 2, 3, Et));		// 行列Etの4行目2列目～3列目を1で埋める (戻り値渡し版)
	constexpr auto Etx4 = fillrow(39, 3, 1, 2, Etx3);	// コンパイル時に指定行を値で埋める
	disp(Etx4);
	ArcsMat<4,1,size_t> vor = {
		2,
		1,
		4,
		3
	};
	disp(vor);
	disp(orderrow(Et, vor));			// 行列Etをvorが昇順になるように行を並び替える (戻り値渡し版のみ)
	disp(orderrow(Et, vor));			// 2回やると元に戻る
	orderrow_and_vec(Et, vor);			// 行列Etをvorが昇順になるようにEtとvorの両方の行を並び替える (引数渡し版のみ)
	disp(Et);
	disp(vor);
	constexpr ArcsMat<4,1,size_t> vorx = {
		2,
		1,
		4,
		3
	};
	constexpr auto Etx5 = orderrow(Etx, vorx);			// コンパイル時に並び替える
	disp(Etx5);
	ArcsMat<4,1> et5;
	sumrow(Etx5, et5);				// 各々の行の総和を計算する (引数渡し版)
	disp(et5);
	disp(sumrow(Etx5));				// 各々の行の総和を計算する (戻り値渡し版)
	constexpr auto etx5 = sumrow(Etx5);					// コンパイル時に総和を計算する
	disp(etx5);
	
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
	disp(F);
	constexpr ArcsMat<4,1> vfx = {
		99,
		88,
		77,
		66
	};
	auto vf = vfx;
	ArcsMat<3,1> vfg;
	getvvector(F, vfg, 3, 2);			// 行列Fの先頭位置(3,2)から縦ベクトルを抽出する (引数渡し版)
	disp(vfg);
	disp(getvvector<3>(F, 2, 3));		// 行列Fの先頭位置(2,3)から縦ベクトルを抽出する (戻り値渡し版)
	constexpr auto vfgx = getvvector<5>(Fx, 2, 5);		// コンパイル時に縦ベクトルを抽出
	disp(vfgx);
	setvvector(F, vf, 3, 5);			// 行列Fの先頭位置(3,5)に縦ベクトルvfを上書きする (引数渡し版)
	disp(F);
	disp(setvvector(vf, 3, 4, F));		// 行列Fの先頭位置(3,4)に縦ベクトルvfを上書きする (戻り値渡し版)
	constexpr auto Fx2 = setvvector(vfx, 3, 4, Fx);		// コンパイル時に縦ベクトルを上書き
	disp(Fx2);
	constexpr ArcsMat<1,4> wfx = { 55, 44, 33, 22 };
	auto wf = wfx;
	ArcsMat<1,3> wfg;
	gethvector(F, wfg, 2, 1);			// 行列Fの先頭位置(2,1)から横ベクトルを抽出する (引数渡し版)
	disp(wfg);
	disp(gethvector<3>(F, 5, 2));		// 行列Fの先頭位置(5,2)から3列の横ベクトルを抽出する (戻り値渡し版)
	constexpr auto wfgx = gethvector<5>(Fx, 6, 2);		// コンパイル時に横ベクトルを抽出
	disp(wfgx);
	sethvector(F, wf, 2, 3);			// 行列Fの先頭位置(2,3)に横ベクトルwfを上書きする (引数渡し版)
	disp(F);
	disp(sethvector(wf, 1, 3, F));		// 行列Fの先頭位置(1,3)に横ベクトルwfを上書きする (戻り値渡し版)
	constexpr auto Fx4 = sethvector(wfx, 2, 3, Fx2);	// コンパイル時に横ベクトルを上書き
	disp(Fx4);
	ArcsMat<4,3> F43;
	getsubmatrix(F, F43, 2, 1);			// 行列Fの左上位置(2,1)から4×3の小行列を抽出する (引数渡し版)
	disp(F43);
	disp(( getsubmatrix<4,3>(F, 2, 2) ));	// 行列Fの左上位置(2,2)から4×3の小行列を抽出する (戻り値渡し版)
	constexpr auto Fx5 = getsubmatrix<2,3>(Fx4, 3, 3);	// コンパイル時に小行列を抽出する
	disp(Fx5);
	F43.Set(
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		0, 1, 2
	);
	setsubmatrix(F, F43, 2, 4);			// 行列Fの左上位置(2,4)に4×3の行列を書き込む (引数渡し版)
	disp(F);
	disp(setsubmatrix(F43, 2, 1, F));	// 行列Fの左上位置(2,1)に4×3の行列を書き込む (戻り値渡し版)
	constexpr ArcsMat<2,3> F23x = {
		100, 200, 300,
		400, 500, 600
	};
	constexpr auto Fx6 = setsubmatrix(F23x, 3, 3, Fx4);	// コンパイル時に小行列を書き込む
	disp(Fx6);

	// 行列間のコピー操作関連の関数
	printf("\n★★★★★★★ 行列間のコピー操作関連の関数\n");
	ArcsMat<10,10> Acpy1;
	copymatrix(Fx, 2,6, 3,4, Acpy1, 5,7);	// 行列Fxから行列Acpy1にコピー (引数渡し版のみ)
	disp(Fx);								// 等価なMATLABコード: Acpy1(5:9, 7:8) = Fx(2:6, 3:4)
	disp(Acpy1);
	
	// シフト関連の関数
	printf("\n★★★★★★★ シフト関連の関数\n");
	F = Fx;
	auto G = ones<7,6>();
	disp(F);
	shiftup(F, G);			// 行列Fを1行分上にシフトして行列Gとして出力 (引数渡し版)
	disp(G);
	shiftup(G, G, 3);		// 行列Gを3行分上にシフト (引数渡し版)
	disp(G);
	disp(shiftup(F));		// 行列Fを1行分上にシフトして出力 (戻り値渡し版)
	disp(shiftup(F, 0));	// 行列Fを0行分上にシフトして出力 (戻り値渡し版)
	disp(shiftup(F, 4));	// 行列Fを4行分上にシフトして出力 (戻り値渡し版)
	constexpr auto Gx1 = shiftup(Fx, 4);				// コンパイル時に上にシフト
	disp(Gx1);
	disp(F);
	shiftdown(F, G);		// 行列Fを1行分下にシフトして行列Gとして出力 (引数渡し版)
	disp(G);
	shiftdown(F, G, 3);		// 行列Fを3行分下にシフトして行列Gとして出力 (引数渡し版)
	disp(G);
	F = Fx;
	disp(F);
	disp(shiftdown(F));		// 行列Fを1行分下にシフトして出力 (戻り値渡し版)
	disp(shiftdown(F, 4));	// 行列Fを4行分下にシフトして出力 (戻り値渡し版)
	constexpr auto Gx2 = shiftdown(Fx, 2);				// コンパイル時に下にシフト
	disp(Gx2);
	F = Fx;
	shiftleft(F, G);			// 行列Fを1列分左にシフトして行列Gとして出力 (引数渡し版)
	disp(G);
	shiftleft(F, G, 3);			// 行列Fを3列分左にシフトして行列Gとして出力 (引数渡し版)
	disp(G);
	disp(shiftleft(F));			// 行列Fを1列分左にシフトして出力 (戻り値渡し版)
	disp(shiftleft(F, 2));		// 行列Fを2列分左にシフトして出力 (戻り値渡し版)
	constexpr auto Gx3 = shiftleft(Fx, 4);				// コンパイル時に左にシフト
	disp(Gx3);
	F = Fx;
	shiftright(F, G);			// 行列Fを1列分右にシフトして行列Gとして出力 (引数渡し版)
	disp(G);
	shiftright(F, G, 3);		// 行列Fを3列分右にシフトして行列Gとして出力 (引数渡し版)
	disp(G);
	disp(shiftright(F));		// 行列Fを1列分右にシフトして出力 (戻り値渡し版)
	disp(shiftright(F, 2));		// 行列Fを2列分右にシフトして出力 (戻り値渡し版)
	constexpr auto Gx4 = shiftright(Fx, 4);				// コンパイル時に右にシフト
	disp(Gx4);
	
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
	disp(Hx1);
	disp(Hx2);
	disp(H12);
	disp(concatv(Hx2, Hx1));	// 行列を縦に連結する (戻り値渡し版)
	constexpr auto Hx12 = concatv(Hx1, Hx2);			// コンパイル時に連結
	disp(Hx12);
	constexpr ArcsMat<4,1> Hx3 = {
		1,
		2,
		3,
		4	
	};
	ArcsMat<4,4> H44;
	concath(Hx2, Hx3, H44);		// 行列を横に連結する (引数渡し版)
	disp(H44);
	disp(concath(Hx2, Hx3));	// 行列を横に連結する (戻り値渡し版)
	constexpr auto Hx44 = concath(Hx2, Hx3);			// コンパイル時に連結
	disp(Hx44);
	constexpr ArcsMat<2,1> Hx4 = {
		5,
		6
	};
	ArcsMat<6,4> H64;
	concat4(Hx1, Hx3, Hx2, Hx4, H64);		// 4つの行列を連結する (引数渡し版)
	disp(H64);
	disp(concat4(Hx1, Hx4, Hx2, Hx3));		// 4つの行列を連結する (戻り値渡し版)
	constexpr auto Hx64 = concat4(Hx1, Hx3, Hx2, Hx4);	// コンパイル時に連結
	disp(Hx64);
	
	// 対角要素操作関連の関数
	printf("\n★★★★★★★ 対角要素操作関連の関数\n");
	constexpr ArcsMat<3,1> jx1 = {1, 3, 5};
	disp(jx1);
	ArcsMat<3,3> Jy1;
	diag(jx1, Jy1);							// 対角要素に縦ベクトルの各要素を持つ行列を生成する (引数渡し版) 
	disp(Jy1);
	disp(diag(jx1));						// 対角要素に縦ベクトルの各要素を持つ行列を生成する (戻り値渡し版) 
	constexpr auto Jx1 = diag(jx1);						// コンパイル時に生成
	disp(Jx1);
	constexpr ArcsMat<3,3> Jx2 = {
		3, 1, 1,
		1, 9, 1,
		1, 1, 2
	};
	disp(Jx2);
	ArcsMat<3,1> jy2;
	getdiag(Jx2, jy2);						// 対角要素を縦ベクトルとして取得する (引数渡し版)
	disp(jy2);
	disp(getdiag(Jx2));						// 対角要素を縦ベクトルとして取得する (戻り値渡し版)
	constexpr auto jx2 = getdiag(Jx2);					// コンパイル時に取得
	disp(jx2);
	constexpr ArcsMat<3,2> Jx3 = {
		3, 5,
		5, 9,
		5, 5
	};
	disp(Jx3);
	disp(getdiag(Jx3));						// 縦長の場合の主対角要素を取得
	disp(getdiag(tp(Jx3)));					// 横長の場合の主対角要素を取得
	disp(B);
	disp(getdiag(B,  1));					// 主対角より1つ上の対角成分を取得
	disp(getdiag(B, -1));					// 主対角より1つ下の対角成分を取得
	disp(Hx64);
	disp(getdiag(Hx64,  2));				// 主対角より2つ上の対角成分を取得
	disp(getdiag(Hx64,  1));				// 主対角より1つ上の対角成分を取得
	disp(getdiag(Hx64, -2));				// 主対角より2つ下の対角成分を取得
	disp(getdiag(Hx64, -3));				// 主対角より3つ下の対角成分を取得
	dispf(~Hx64, "% 4.0f");
	disp(getdiag(~Hx64,  3));				// 主対角より3つ上の対角成分を取得
	disp(getdiag(~Hx64,  2));				// 主対角より2つ上の対角成分を取得
	disp(getdiag(~Hx64, -1));				// 主対角より1つ下の対角成分を取得
	disp(getdiag(~Hx64, -2));				// 主対角より2つ下の対角成分を取得
	printf("trace(Jx2) = %f\n", trace(Jx2));// トレースを計算する (戻り値渡し版のみ)
	constexpr double trJx2 = trace(Jx2);				// コンパイル時にトレースを計算
	printf("trJx2 = %f\n", trJx2);
	printf("multdiag(Jx2) = %f\n", multdiag(Jx2));		// 対角要素の総積を計算する (戻り値渡し版のみ)
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
	disp(k1);
	disp(Fx);
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
	using namespace std::literals::complex_literals;	// 虚数単位リテラル「i」の使用
	ArcsMat<1,3,std::complex<double>> k1comp = {-2.0 + 2.0i, 4.0 + 1.0i, -1.0 - 3.0i};
	disp(k1comp);
	auto k1cmax = max(k1comp);					// 複素数の最大値
	printf("max(k1comp) = %f + %fj\n", std::real(k1cmax), std::imag(k1cmax));
	auto k1cmin = min(k1comp);					// 複素数の最小値
	printf("min(k1comp) = %f + %fj\n", std::real(k1cmin), std::imag(k1cmin));

	// 要素ごとの数学関数
	printf("\n★★★★★★★ 要素ごとの数学関数\n");
	constexpr ArcsMat<3,4> Ax1 = {
		       1.0,        2.0,        3.0,      4.0,
		    M_PI_4,     M_PI_2, 3.0*M_PI_4,	    M_PI,
		5.0*M_PI_4, 3.0*M_PI_2, 7.0*M_PI_4, 2.0*M_PI 
	};
	dispf(Ax1, "%8.3f");
	printf("sum(Ax1) = %f\n", sum(Ax1));	// 行列要素の総和 (戻り値渡し版のみ)
	ArcsMat<3,4> Y1;
	exp(Ax1, Y1);							// 行列要素の指数関数 (引数渡し版)
	dispf(exp(Ax1), "%8.3f");				// 行列要素の指数関数 (戻り値渡し版)
	constexpr auto Yx1 = exp(Ax1);						// コンパイル時に行列要素の指数関数を計算
	dispf(Yx1, "%8.3f");
	ArcsMat<3,4> Y2;
	log(Y1, Y2);							// 行列要素の対数関数(底e版) (引数渡し版)
	dispf(log(Y1), "%8.3f");				// 行列要素の対数関数(底e版) (戻り値渡し版)
	constexpr auto Yx2 = log(Yx1);						// コンパイル時に行列要素の対数関数(底e版)を計算
	dispf(Yx2, "%8.3f");
	log10(Ax1, Y1);							// 行列要素の対数関数(底10版) (引数渡し版)
	dispf(log10(Ax1), "%8.3f");				// 行列要素の対数関数(底10版) (戻り値渡し版)
	constexpr auto Yx3 = log10(Ax1);					// コンパイル時に行列要素の対数関数(底10版)を計算
	dispf(Yx3, "%8.3f");
	sin(Ax1, Y1);							// 行列要素の正弦関数 (引数渡し版)
	dispf(sin(Ax1), "%8.3f");				// 行列要素の正弦関数 (戻り値渡し版)
	constexpr auto Yx4 = sin(Ax1);						// コンパイル時に行列要素の正弦関数を計算
	dispf(Yx4, "%8.3f");
	cos(Ax1, Y1);							// 行列要素の余弦関数 (引数渡し版)
	dispf(cos(Ax1), "%8.3f");				// 行列要素の余弦関数 (戻り値渡し版)
	constexpr auto Yx5 = cos(Ax1);						// コンパイル時に行列要素の余弦関数を計算
	dispf(Yx5, "%8.3f");
	tan(Ax1, Y1);							// 行列要素の正接関数 (引数渡し版)
	dispf(tan(Ax1), "% 12.3e");				// 行列要素の正接関数 (戻り値渡し版)
	constexpr auto Yx6 = tan(Ax1);						// コンパイル時に行列要素の正接関数を計算
	dispf(Yx6, "% 12.3e");
	sqrt(Ax1, Y1);							// 行列要素の平方根 (引数渡し版)
	dispf(sqrt(Ax1), "% 8.3f");				// 行列要素の平方根 (戻り値渡し版)
	constexpr auto Yx7 = sqrt(Ax1);						// コンパイル時に行列要素の平方根を計算
	dispf(Yx7, "% 8.3f");
	sqrt(Ax1, Y1);							// 行列要素の符号関数 (引数渡し版)
	dispf(sqrt(Ax1), "% 8.3f");				// 行列要素の符号関数 (戻り値渡し版)
	constexpr auto Yx8 = sign(Ax1);						// コンパイル時に行列要素の符号関数を計算
	dispf(Yx8, "% 8.3f");
	abs(Yx5, Y1);							// 行列要素の絶対値 (引数渡し版)
	dispf(abs(Yx5), "% 12.3f");				// 行列要素の絶対値 (戻り値渡し版)
	constexpr auto Yx9 = abs(Yx5);						// コンパイル時に行列要素の絶対値を計算
	dispf(Yx9, "% 12.3f");

	// 複素数関連の関数
	printf("\n★★★★★★★ 複素数関連の関数\n");
	constexpr ArcsMat<3,5,std::complex<double>> Acmpx1 = {
		-2, -1,  0,  1,  2,
		-4, -3,  0,  3,  4,
		-1, -2, -1, -2, -1
	};
	disp(Acmpx1);
	ArcsMat<3,5,std::complex<double>> Y9;
	sqrt(Acmpx1, Y9);						// 負の行列要素の平方根 (引数渡し版)
	dispf(sqrt(Acmpx1), "% 6.3f");			// 負の行列要素の平方根 (戻り値渡し版)
	//constexpr auto Yx9 = sqrt(Acmpx1);			// コンパイル時に行列要素の平方根を計算 驚異の非対応エラー!!(C++20以上にしないと解決不可)
	//dispf(Yx9, "%8.3f");
	constexpr ArcsMat<2,3,std::complex<double>> Acmpx2 = {
		std::complex( 1.0, 1.0), std::complex( 3.0, 4.0), std::complex( 3.0,-4.0),
		std::complex(-1.0, 1.0), std::complex(-3.0, 4.0), std::complex(-3.0,-4.0)
	};
	disp(Acmpx2);
	ArcsMat<2,3,double> Y10;
	abs(Acmpx2, Y10);						// 複素数行列要素の絶対値 (引数渡し版)
	disp(abs(Acmpx1));						// 複素数行列要素の絶対値 (戻り値渡し版)
	arg(Acmpx2, Y10);							// 複素数行列要素の偏角 (引数渡し版)   [rad]
	dispf(arg(Acmpx2)*180.0/M_PI, "% 6.1f");	// 複素数行列要素の偏角 (戻り値渡し版) [deg]
	real(Acmpx2, Y10);						// 複素数行列要素の実数部 (引数渡し版)
	disp(real(Acmpx2));						// 複素数行列要素の実数部 (戻り値渡し版)
	imag(Acmpx2, Y10);						// 複素数行列要素の実数部 (引数渡し版)
	disp(imag(Acmpx2));						// 複素数行列要素の実数部 (戻り値渡し版)
	ArcsMat<2,3,std::complex<double>> Y11;
	conj(Acmpx2, Y11);						// 複素数行列要素の複素共役 (引数渡し版)
	disp(conj(Acmpx2));						// 複素数行列要素の複素共役 (戻り値渡し版)
	ArcsMat<3,2,std::complex<double>> Y12;
	Y9.FillAll(39);							// 複素数行列の実数部に値を埋める
	disp(Y9);
	Y9.FillAll( 3.9 + 2.4i );				// 複素数行列に値を埋める
	disp(Y9);
	ArcsMat<2,3,std::complex<double>> Acmp1 = Pi;		//「実数行列 → 複素数行列」のコピーコンストラクタ
	disp(Acmp1);
	
	// 転置行列関連の関数と演算子
	printf("\n★★★★★★★ 転置行列関連の関数と演算子\n");
	ArcsMat<2,3> Dt;
	disp(D);
	tp(D, Dt);							// Dの転置 (引数渡し版)
	disp(tp(D));						// Dの転置 (戻り値渡し版)
	constexpr auto Atx = tp(Ax);		// コンパイル時に転置行列を生成
	disp(Ax);
	disp(Atx);
	//tp(A, Dt);	// ←サイズエラー！(アンコメントするとAssertion Failed)
	//A = tp(Dt);	// ←サイズエラー！(アンコメントするとAssertion Failed)
	Htp(Acmpx2, Y12);					// エルミート転置 (引数渡し版)
	disp(Htp(Acmpx2));					// エルミート転置 (戻り値渡し版)
	disp(~D);							// 実数行列の転置 (演算子版)
	disp(~Acmpx2);						// 複素数行列のエルミート転置 (演算子版は自動的に共役転置)
	constexpr auto jxjtx = jx1*~jx1;	// コンパイル時に転置して乗算 (演算子版)
	constexpr auto jxtjx = (~jx1)*jx1;	// コンパイル時に転置して乗算 (演算子版)
	disp(jx1);
	dispf(jxjtx, "% 3.0f");
	dispf(jxtjx, "% 3.0f");

	// ノルム関連の関数(行列版)
	printf("\n★★★★★★★ ノルム関連の関数(行列版)\n");
	dispf(Ax1, "% 7.4f");
	printf("norm<euc>(Ax1) = % 7.4f\n", norm<NormType::AMT_L2>(Ax1));	// ユークリッドL2ノルムを計算する (戻り値渡し版のみ)
	double enormAx1 = norm<NormType::AMT_L2>(Ax1);						// コンパイル時にユークリッドL2ノルムを計算
	printf("norm<euc>(Ax1) = % 7.4f\n", enormAx1);
	printf("norm<man>(Ax1) = % 7.4f\n", norm<NormType::AMT_L1>(Ax1));	// 絶対値L1ノルムを計算する (戻り値渡し版のみ)
	constexpr double mnormAx1 = norm<NormType::AMT_L1>(Ax1);			// コンパイル時に絶対値L1ノルムを計算
	printf("norm<man>(Ax1) = % 7.4f\n", mnormAx1);
	printf("norm<inf>(Ax1) = % 7.4f\n", norm<NormType::AMT_LINF>(Ax1));	// 無限大L∞ノルムを計算する (戻り値渡し版のみ)
	constexpr double inormAx1 = norm<NormType::AMT_LINF>(Ax1);			// コンパイル時に無限大L∞ノルムを計算
	printf("norm<inf>(Ax1) = % 7.4f\n", inormAx1);
	dispf(Acmpx2, "% 7.4f");
	printf("norm<euc>(Acmpx2) = % 7.4f\n", norm<NormType::AMT_L2>(Acmpx2));		// 複素数ユークリッドL2ノルムを計算する (戻り値渡し版のみ)
	printf("norm<man>(Acmpx2) = % 7.4f\n", norm<NormType::AMT_L1>(Acmpx2));		// 複素数絶対値L1ノルムを計算する (戻り値渡し版のみ)
	printf("norm<inf>(Acmpx2) = % 7.4f\n", norm<NormType::AMT_LINF>(Acmpx2));	// 複素数無限大L∞ノルムを計算する (戻り値渡し版のみ)
	
	// ノルム関連の関数(ベクトル版)
	printf("\n★★★★★★★ ノルム関連の関数(ベクトル版)\n");
	dispf(k1, "% 7.4f");
	printf("norm<euc>(k1) = % 7.4f\n", norm<NormType::AMT_L2>(k1));		// ユークリッドL2ノルムを計算する (戻り値渡し版のみ)
	constexpr double enormk1 = norm<NormType::AMT_L2>(k1);				// コンパイル時にユークリッドL2ノルムを計算
	printf("norm<euc>(k1) = % 7.4f\n", enormk1);
	printf("norm<man>(k1) = % 7.4f\n", norm<NormType::AMT_L1>(k1));		// 絶対値L1ノルムを計算する (戻り値渡し版のみ)
	constexpr double mnormk1 = norm<NormType::AMT_L1>(k1);				// コンパイル時に絶対値L1ノルムを計算
	printf("norm<man>(k1) = % 7.4f\n", mnormk1);
	printf("norm<inf>(k1) = % 7.4f\n", norm<NormType::AMT_LINF>(k1));	// 無限大L∞ノルムを計算する (戻り値渡し版のみ)
	constexpr double inormk1 = norm<NormType::AMT_LINF>(k1);			// コンパイル時に無限大L∞ノルムを計算
	printf("norm<inf>(k1) = % 7.4f\n\n", inormk1);
	auto k1cmpx = sqrt(static_cast<ArcsMat<7,1,std::complex<double>>>(-k1));
	dispf(k1cmpx, "% 7.4f");
	printf("norm<euc>(k1cmpx) = % 7.4f\n", norm<NormType::AMT_L2>(k1cmpx));		// 複素数ユークリッドL2ノルムを計算する (戻り値渡し版のみ)
	printf("norm<man>(k1cmpx) = % 7.4f\n", norm<NormType::AMT_L1>(k1cmpx));		// 複素数絶対値L1ノルムを計算する (戻り値渡し版のみ)
	printf("norm<inf>(k1cmpx) = % 7.4f\n", norm<NormType::AMT_LINF>(k1cmpx));	// 複素数無限大L∞ノルムを計算する (戻り値渡し版のみ)

	// 三角行列操作関連の関数
	printf("\n★★★★★★★ 三角行列操作関連の関数\n");
	ArcsMat<7,6> Fx7;
	gettriup(Fx, Fx7);						// 上三角行列を切り出す (引数渡し版)
	dispf(Fx7, "%3.0f");
	dispf(gettriup(Fx, 3), "%3.0f");		// 上三角行列を切り出す (戻り値渡し版)
	constexpr auto Fxtri = gettriup(Fx, 2);				// コンパイル時に上三角行列を切り出す 
	dispf(Fxtri, "%3.0f");
	Fx7.FillAllZero();
	gettrilo(Fx, Fx7);						// 下三角行列を切り出す (引数渡し版)
	dispf(Fx7, "%3.0f");
	dispf(gettrilo(Fx, 3), "%3.0f");		// 下三角行列を切り出す (戻り値渡し版)
	constexpr auto Fxtri2 = gettrilo(Fx, 2);			// コンパイル時に下三角行列を切り出す 
	dispf(Fxtri2, "%3.0f");

	// LU分解関連の関数
	printf("\n★★★★★★★ LU分解関連の関数\n");
	ArcsMat<3,3> L, U, P;
	A.Set(
		10, -7,  0,
		-3,  2,  6,
		 5, -1,  5
	);
	dispf(A, "% 3.0f");
	LUP(A, L, U, P);				// LU分解の結果と置換行列Pを計算 (引数渡し版)
	dispf(L, "% 6.2f");
	dispf(U, "% 6.2f");
	dispf(P, "% 3.0f");
	dispf(~P*L*U, "% 3.0f");		// もとに戻るかチェック
	std::tie(L, U, P) = LUP(A);		// LU分解の結果と置換行列Pを計算 (タプル返し版)
	dispf(L, "% 6.2f");
	dispf(U, "% 6.2f");
	dispf(P, "% 3.0f");
	LU(A, L, U);					// LU分解の結果のみを計算 (引数渡し版)
	dispf(L, "% 6.2f");
	dispf(U, "% 6.2f");
	dispf(L*U, "% 3.0f");			// もとに戻るかチェック
	std::tie(L, U) = LU(A);			// LU分解の結果のみを計算 (タプル返し版)
	dispf(L, "% 6.2f");
	dispf(U, "% 6.2f");
	dispf(L*U, "% 3.0f");			// もとに戻るかチェック
	constexpr auto LxUx = LU(Ax);						// コンパイル時にLU分解を計算
	constexpr auto Lx = std::get<0>(LxUx);				// コンパイル時に計算した下三角を抽出
	constexpr auto Ux = std::get<1>(LxUx);				// コンパイル時に計算した上三角を抽出
	dispf(Ax, "% 3.0f");
	dispf(Lx, "% 6.3f");
	dispf(Ux, "% 6.3f");
	dispf(Lx*Ux, "%3.0f");			// もとに戻るかチェック
	ArcsMat<3,3,std::complex<double>> Acomp1 = {
		4.0 + 6.0i,  1.0 - 3.0i,  5.0 + 2.0i,
		8.0 - 5.0i, -7.0 - 6.0i,  7.0 - 1.0i,
		9.0 + 9.0i, -7.0 - 5.0i, -5.0 - 3.0i
	};
	dispf(Acomp1, "% 3.0f");
	auto [Lcomp, Ucomp, Pcomp] = LUP(Acomp1);	// 複素数LU分解と置換行列の計算
	dispf(Lcomp, "% 8.3f");
	dispf(Ucomp, "% 8.3f");
	dispf(Pcomp, "% 8.3f");
	dispf(~Pcomp*Lcomp*Ucomp, "% 5.2f");		// もとに戻るかチェック

	// 行列式det関連の関数
	printf("\n★★★★★★★ 行列式det関連の関数\n");
	ArcsMat<3,3> Adet1 = {
		 1, -2,  4,
    	-5,  2,  0,
     	 1,  0,  3
	};
	printf("det(Adet1) = % 8.4f\n", det(Adet1));// 行列式の計算1
	printf("det(A) = % 8.4f\n", det(A));		// 行列式の計算2
	printf("det(Ax) = % 8.4f\n", det(Ax));		// 行列式の計算3
	constexpr double detAx = det(Ax);					// コンパイル時に行列式を計算
	printf("det(Ax) = % 8.4f\n", detAx);				// コンパイル時に計算した行列式を表示
	std::complex detAcomp1 = det(Acomp1);		// 複素数の行列式の計算
	printf("det(Acomp1) = % 8.4f + % 8.4fi\n", std::real(detAcomp1), std::imag(detAcomp1));

	// Householder行列関連の関数
	printf("\n★★★★★★★ Householder行列関連の関数\n");
	constexpr ArcsMat<3,1> vhld1 = {9, 3, 0};
	dispf(vhld1, "% 8.4f");
	ArcsMat<3,3> Hhld1;
	Householder(vhld1, Hhld1);		// ハウスホルダー行列を生成 (引数渡し版)
	Hhld1 = Householder(vhld1);		// ハウスホルダー行列を生成 (戻り値返し版)
	dispf(Hhld1, "% 8.4f");
	Householder(vhld1, Hhld1, 2);	// 次元を2にする場合 (引数渡し版)
	Hhld1 = Householder(vhld1, 2);	// 次元を2にする場合 (戻り値返し版)
	dispf(Hhld1, "% 8.4f");
	constexpr auto Hhld1x = Householder(vhld1);			// コンパイル時にハウスホルダー行列を生成
	dispf(Hhld1x, "% 8.4f");
	
	// QR分解関連の関数
	printf("\n★★★★★★★ QR分解関連の関数\n");
	constexpr ArcsMat<5,5> Aqr1 = {
		17, 24,  1,  8, 15,
		23,  5,  7, 14, 16,
		 4,  6, 13, 20, 22,
		10, 12, 19, 21,  3,
		11, 18, 25,  2,  9
	};
	ArcsMat<5,5> Qqr1, Rqr1;
	QR(Aqr1, Qqr1, Rqr1);				// QR分解を計算 (引数渡し版)
	dispf(Aqr1, "% 8.4f");
	dispf(Qqr1, "% 8.4f");
	dispf(Rqr1, "% 8.4f");
	dispf(Qqr1*~Qqr1, "% 8.4f");		// Qが直交行列かチェック
	dispf(Qqr1*Rqr1, "% 8.4f");			// 元に戻るかチェック
	printf("||Aqr1 - Qqr1*Rqr1|| = %e\n\n", norm(Aqr1 - Qqr1*Rqr1));	// ユークリッドL2ノルムでチェック
	constexpr ArcsMat<3,4> Aqr2 = {
		12, -51,   4, 39,
		 6, 167, -68, 22,
		-4,  24, -41, 11
	};
	auto [Qqr2, Rqr2] = QR(Aqr2);		// QR分解を計算 (タプル返し版)
	dispf(Aqr2, "% 5.0f");
	dispf(Qqr2, "% 8.4f");
	dispf(Rqr2, "% 8.4f");
	dispf(Qqr2*~Qqr2, "% 8.4f");		// Qが直交行列かチェック
	dispf(Qqr2*Rqr2, "% 8.4f");			// 元に戻るかチェック
	constexpr auto QRx2 = QR(~Aqr2);					// コンパイル時にQR分解を計算
	constexpr auto Qx2 = std::get<0>(QRx2);				// コンパイル時に計算したユニタリ行列を抽出
	constexpr auto Rx2 = std::get<1>(QRx2);				// コンパイル時に計算したユニタリ行列を抽出
	dispf(~Aqr2, "% 5.0f");
	dispf(Qx2, "% 8.4f");
	dispf(Rx2, "% 8.4f");
	dispf(Qx2*~Qx2, "% 8.4f");			// Qが直交行列かチェック
	dispf(Qx2*Rx2, "% 8.4f");			// 元に戻るかチェック
	auto [Qqr3, Rqr3] = QR(Acomp1);		// 複素数QR分解を計算 正方行列 (タプル返し版)
	dispf(Acomp1, "% 2.0f");
	dispf(Qqr3, "% 8.4f");
	dispf(Rqr3, "% 8.4f");
	dispf(Qqr3*~Qqr3, "% 3.1f");		// Qが直交行列かチェック
	dispf(Qqr3*Rqr3, "% 3.1f");			// 元に戻るかチェック
	auto [Qqr4, Rqr4] = QR(Acmpx2);		// 複素数QR分解を計算 横長行列 (タプル返し版)
	dispf(Acmpx2, "% 2.0f");
	dispf(Qqr4, "% 8.4f");
	dispf(Rqr4, "% 8.4f");
	dispf(Qqr4*~Qqr4, "% 3.1f");		// Qが直交行列かチェック
	dispf(Qqr4*Rqr4, "% 3.1f");			// 元に戻るかチェック
	auto [Qqr5, Rqr5] = QR(~Acmpx2);	// 複素数QR分解を計算 縦長行列 (タプル返し版)
	dispf(~Acmpx2, "% 2.0f");
	dispf(Qqr5, "% 8.4f");
	dispf(Rqr5, "% 8.4f");
	dispf(Qqr5*~Qqr5, "% 3.1f");		// Qが直交行列かチェック
	dispf(Qqr5*Rqr5, "% 3.1f");			// 元に戻るかチェック

	// SVD特異値分解関連の関数
	printf("\n★★★★★★★ SVD特異値分解関連の関数\n");
	constexpr ArcsMat<4,2> As1 = {
		1, 2,
		3, 4,
		5, 6,
		7, 8
	};
	ArcsMat<4,4> Us1;
	ArcsMat<4,2> Ss1;
	ArcsMat<2,2> Vs1;
	SVD(As1, Us1, Ss1, Vs1);			// 特異値分解 (引数渡し版)：縦長行列の場合
	dispf(As1, "% 8.4f");
	dispf(Us1, "% 8.4f");
	dispf(Ss1, "% 8.4f");
	dispf(Vs1, "% 8.4f");
	dispf(Us1*Ss1*~Vs1, "% 8.4f");		// 元に戻るかチェック
	constexpr auto USVs2 = SVD(~As1);					// コンパイル時に特異値分解を計算 (タプル返し版)：横長行列の場合
	constexpr auto Us2 = std::get<0>(USVs2);			// コンパイル時に計算したU行列を抽出
	constexpr auto Ss2 = std::get<1>(USVs2);			// コンパイル時に計算したΣ行列を抽出
	constexpr auto Vs2 = std::get<2>(USVs2);			// コンパイル時に計算したV行列を抽出
	dispf(Us2, "% 8.4f");
	dispf(Ss2, "% 8.4f");
	dispf(Vs2, "% 8.4f");
	dispf(Us2*Ss2*~Vs2, "% 8.4f");		// 元に戻るかチェック
	ArcsMat<3,3> As2 = {
		 2,  0,  2,
		 0,  1,  0,
		 0,  0,  0
	};
	ArcsMat<3,3> Us3, Ss3, Vs3;
	std::tie(Us3, Ss3, Vs3) = SVD(As2);	// 特異値分解 (タプル返し版)：ランク落ちの場合
	dispf(Us3, "% 8.4f");
	dispf(Ss3, "% 8.4f");
	dispf(Vs3, "% 8.4f");
	dispf(Us3*Ss3*~Vs3, "% 8.3f");		// 元に戻るかチェック
	As2.Set(
		 1,  1,  3,
		-5,  6, -3, 
		 7, -2,  9
	);
	std::tie(Us3, Ss3, Vs3) = SVD(As2);	// 特異値分解 (タプル返し版)：符号修正が必要な場合
	dispf(As2, "% 8.4f");
	dispf(Us3, "% 8.4f");
	dispf(Ss3, "% 8.4f");
	dispf(Vs3, "% 8.4f");
	dispf(Us3*Ss3*~Vs3, "% 8.4f");		// 元に戻るかチェック
	ArcsMat<3,3,std::complex<double>> Us4, Ss4, Vs4;
	std::tie(Us4, Ss4, Vs4) = SVD(Acomp1);	// 複素数特異値分解を計算 (タプル返し版)
	dispf(Acomp1, "% 8.4f");
	dispf(Us4, "% 8.4f");
	dispf(Ss4, "% 8.4f");
	dispf(Vs4, "% 8.4f");
	dispf(Us4*Ss4*~Vs4, "% 8.4f");		// 元に戻るかチェック
	auto [Us5, Ss5, Vs5] = SVD(Acmpx2);	// 複素数特異値分解を計算 横長行列 (タプル返し版)
	dispf(Acmpx2, "% 2.0f");
	dispf(Us5, "% 8.4f");
	dispf(Ss5, "% 8.4f");
	dispf(Vs5, "% 8.4f");
	dispf(Us5*Ss5*~Vs5, "% 8.4f");		// 元に戻るかチェック
	ArcsMat<2,5,std::complex<double>> Acomp2 = {
		7.0000 + 8.0000i,   5.0000 + 2.0000i,   7.0000 + 5.0000i,   1.0000 + 6.0000i,   1.0000 + 5.0000i,
		5.0000 + 7.0000i,   1.0000 + 6.0000i,   1.0000 + 9.0000i,   5.0000 + 8.0000i,   8.0000 + 4.0000i
	};
	auto [Us6, Ss6, Vs6] = SVD(Acomp2);	// 複素数特異値分解を計算 長い横長行列 (タプル返し版)
	dispf(Acomp2, "% 8.4f");
	dispf(Us6, "% 8.4f");
	dispf(Ss6, "% 8.4f");
	dispf(Vs6, "% 8.4f");
	dispf(Us6*Ss6*~Vs6, "% 8.4f");		// 元に戻るかチェック

	// 階数関連の関数
	printf("\n★★★★★★★ 階数関連の関数\n");
	constexpr ArcsMat<3,3> Arank1 = {
		 3, 2,  4,
		-1, 1,  2,
		 9, 5, 10
	};
	dispf(Arank1, "% 8.4f");	
	printf("rank(Arank1) = %zu\n", rank(Arank1));	// 階数を計算(戻り値返し版のみ)
	constexpr size_t rankA1 = rank(Arank1);				// コンパイル時にランクを計算
	printf("rank(Arank1) = %zu\n\n", rankA1);
	ArcsMat<4,4> Arank2 = {
		10.0,  0.0,  0.0, 0.0,
		 0.0, 25.0,  0.0, 0.0,
		 0.0,  0.0, 34.0, 0.0,
		 0.0,  0.0,  0.0, 1e-15
	};
	dispf(Arank2, "% 8.4f");	
	printf("rank(Arank2) = %zu\n", rank(Arank2));		// 要素の1つが 1e-15 で小さすぎてランク落ちと見なされる
	printf("rank(Arank2) = %zu\n", rank(Arank2, 1e-16));// そこで、許容誤差を 1e-16 に設定するとフルランクと見なされる

	// コレスキー分解関連の関数
	// 注意：コレスキー分解は対称正定値でないと失敗する(理論的に)が、関数内部での判定は行っていない。
	printf("\n★★★★★★★ コレスキー分解関連の関数\n");
	constexpr ArcsMat<3,3> Achol1 = {
		 7, -2,  9,
		-2,  6, -3,
		 9, -3,  3
	};
	ArcsMat<3,3> Lldl1, Dldl1;
	LDL(Achol1, Lldl1, Dldl1);				// LDL分解を計算 (引数渡し版)
	std::tie(Lldl1, Dldl1) = LDL(Achol1);	// LDL分解を計算 (タプル返し版)
	dispf(Achol1, "% 8.4f");
	dispf(Lldl1, "% 8.4f");
	dispf(Dldl1, "% 8.4f");
	dispf(Lldl1*Dldl1*~Lldl1, "% 8.4f");	// 元に戻るかチェック
	constexpr auto LDldl1x = LDL(Achol1);				// コンパイル時にLDL分解を計算
	constexpr auto Lldl1x = std::get<0>(LDldl1x);		// コンパイル時に計算したL行列を抽出
	constexpr auto Dldl1x = std::get<1>(LDldl1x);		// コンパイル時に計算したL行列を抽出
	dispf(Lldl1x, "% 8.4f");
	dispf(Dldl1x, "% 8.4f");
	ArcsMat<3,3,std::complex<double>> Acomp3 = {
		4.0 + 6.0i,  1.0 - 3.0i,  5.0 + 2.0i,
		1.0 - 3.0i, -7.0 - 6.0i,  7.0 - 1.0i,
		5.0 + 2.0i,  7.0 - 1.0i, -5.0 - 3.0i
	};
	auto [Lldl2, Dldl2] = LDL(Acomp3);		// 複素数LDL分解を計算(MATLABでは非対応で失敗する可能性アリ) (タプル返し版)
	dispf(Acomp3, "% 8.4f");
	dispf(Lldl2, "% 8.4f");
	dispf(Dldl2, "% 8.4f");
	dispf(Lldl2*Dldl2*tp(Lldl2), "% 8.4f");	// 共役転置ではなく普通の転置で元に戻る
	ArcsMat<3,3> Rchol1;
	Cholesky(Achol1, Rchol1);				// 修正コレスキー分解を計算 (引数渡し版)
	dispf(Achol1, "% 8.4f");				// ←対称正定値行列ではないので、
	dispf(Rchol1, "% 8.4f");				// ←NaNが出現
	constexpr ArcsMat<3,3> Achol2 = {
		 2.3370,   -0.0127,   -0.6299,
		-0.0127,    2.3370,   -0.0127,
		-0.6299,   -0.0127,    2.3370
	};
	Cholesky(Achol2, Rchol1);				// 修正コレスキー分解を計算 (引数渡し版)
	Rchol1 = Cholesky(Achol2);				// 修正コレスキー分解を計算 (戻り値返し版)
	dispf(Achol2, "% 8.4f");				// ←対称正定値行列なので、
	dispf(Rchol1, "% 8.4f");				// ←ちゃんと分解できる
	dispf(~Rchol1*Rchol1, "% 8.4f");		// 元に戻るかチェック
	constexpr auto Rchol1x = Cholesky(Achol2);			// コンパイル時に修正コレスキー分解を計算
	dispf(Rchol1x, "% 8.4f");

	// 線形方程式関連の関数
	printf("\n★★★★★★★ 線形方程式関連の関数\n");
	ArcsMat<3,3> Aslv1 = {
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	};
	ArcsMat<3,1> bslv1 = {
		9,
		5,
		7
	};
	ArcsMat<3,1> xslv1;
	linsolve(Aslv1, bslv1, xslv1);	// Ax = bの形の線形方程式をxについて解く関数(引数渡し版)
	dispf(Aslv1, "% 8.4f");			// Aが正方行列で、
	dispf(bslv1, "% 8.4f");			// bがベクトルで、
	dispf(xslv1, "% 8.4f");			// xもベクトルの場合
	ArcsMat<3,3> Bslv2 = {
		9,  2,  7,
		5,  6, -7,
		7, -3,  5
	};
	auto Xslv2 = linsolve(Aslv1, Bslv2);	// AX = Bの形の線形方程式をXについて解く関数(戻り値返し版)
	dispf(Bslv2, "% 8.4f");			// Aが正方行列で、Bが行列で、
	dispf(Xslv2, "% 8.4f");			// Xも行列の場合
	constexpr ArcsMat<4,3> Aslv3 = {
		1,  1,  1,
		2,  3, -2,
		3, -1,  1,
		7,  5,  4
	};
	constexpr ArcsMat<4,1> bslv3 = {
		9,
		5,
		7,
		3
	};
	constexpr auto xslv3 = linsolve(Aslv3, bslv3);	// コンパイル時にAx = bの形の線形方程式をXについて解く
	dispf(Aslv3, "% 8.4f");			// Aが非正方縦長行列で、
	dispf(bslv3, "% 8.4f");			// bがベクトルで、
	dispf(xslv3, "% 8.4f");			// xもベクトルの場合
	ArcsMat<4,3> Bslv4 = {
		9,  2,  7,
		5,  6, -7,
		7, -3,  5,
		3,  1,  2 
	};
	ArcsMat<3,3> Xslv4;
	linsolve(Aslv3, Bslv4, Xslv4);	// AX = Bの形の線形方程式をXについて解く関数(引数渡し版)
	dispf(Bslv4, "% 8.4f");			// Aが非正方縦長行列で、Bが行列で、
	dispf(Xslv4, "% 8.4f");			// Xも行列の場合
	ArcsMat<4,1> xslv5;
	xslv5 = linsolve(~Aslv3, bslv1);// Ax = bの形の線形方程式をxについて解く関数(引数渡し版)
	dispf(xslv5, "% 8.3f");			// Aが非正方横長行列で、bがベクトルで、xもベクトルの場合
	dispf(~Aslv3*xslv5, "% 8.4f");	// この関数はMATLABとは異なる解を出力するが、ただしもちろん Ax = b は成立
	ArcsMat<4,4> Xslv6;				
	Xslv6 = linsolve(~Aslv3, ~Bslv4);	// AX = Bの形の線形方程式をXについて解く関数(引数渡し版)
	dispf(Xslv6, "% 8.3f");				// Aが非正方横長行列で、Bが行列で、Xも行列の場合
	dispf(~Aslv3*Xslv6, "% 8.4f");		// この関数はMATLABとは異なる解を出力するが、ただしもちろん Ax = b は成立

	// 線形最小二乗法関連の関数(linsolveの応用)
	printf("\n★★★★★★★ 線形最小二乗法関連の関数(linsolveの応用)\n");
	const ArcsMat<6,1> tlsq1 = { 0.0, 0.3, 0.8, 1.1, 1.6, 2.3 };			// 説明変数
	const ArcsMat<6,1> ylsq1 = { 0.82, 0.72, 0.63, 0.60, 0.55, 0.50 };		// 観測データ
	const ArcsMat<6,2> Alsq1 = concath(ArcsMat<6,1>::ones(), exp(-tlsq1));	// y(t) = c1 + c2*exp(-t) の行列データ
	dispf(tlsq1, "% 8.4f");
	dispf(ylsq1, "% 8.4f");
	dispf(Alsq1, "% 8.4f");
	auto xlsq1 = linsolve(Alsq1, ylsq1);	// 線形最小二乗法(優決定問題をlinsolveで解くだけ)
	dispf(xlsq1, "% 8.4f");					// x = [ c1, c2 ]^T のフィッティング係数
	printf("||Ax - y|| = % 8.4f\n", norm(Alsq1*xlsq1 - ylsq1));

	// 逆行列関連の関数
	printf("\n★★★★★★★ 逆行列関連の関数\n");
	constexpr ArcsMat<3,3> Ainv1 = {
		 1,  0,  2,
		-1,  5,  0,
		 0,  3, -9
	};
	ArcsMat<3,3> Yinv1;
	inv(Ainv1, Yinv1);			// 逆行列を計算 (引数渡し版)
	Yinv1 = inv(Ainv1);			// 逆行列を計算 (戻り値返し版)
	dispf(Yinv1, "% 8.4f");
	constexpr auto Yinv2 = inv(Ainv1);			// コンパイル時に逆行列を計算
	dispf(Yinv2, "% 8.4f");
	Yinv1 = Ainv1^(-1);			// -1乗と逆行列は等価
	dispf(Yinv1, "% 8.4f");
	auto Yinv3 = inv(Acomp1);	// 複素数逆行列を計算
	dispf(Yinv3, "% 8.4f");
	ArcsMat<3,4> Yinv4;
	pinv(Aslv3, Yinv4);			// Moore-Penroseの擬似逆行列を計算(縦長→横長) (引数渡し版)
	Yinv4 = pinv(Aslv3);		// Moore-Penroseの擬似逆行列を計算(縦長→横長) (戻り値返し版)
	dispf(Yinv4, "% 8.4f");
	constexpr auto Yinv4x = pinv(Aslv3);		// コンパイル時にMoore-Penroseの擬似逆行列を計算
	dispf(Yinv4x, "% 8.4f");
	ArcsMat<4,3> Yinv5;
	pinv(~Aslv3, Yinv5);		// Moore-Penroseの擬似逆行列を計算(横長→縦長) (引数渡し版)
	Yinv5 = pinv(~Aslv3);		// Moore-Penroseの擬似逆行列を計算(横長→縦長) (戻り値返し版)
	dispf(Yinv5, "% 8.4f");
	constexpr auto Yinv5x = pinv(~Aslv3);		// コンパイル時にMoore-Penroseの擬似逆行列を計算
	dispf(Yinv5x, "% 8.4f");

	// Hessenberg分解関連の関数
	printf("\n★★★★★★★ Hessenberg分解関連の関数\n");
	constexpr ArcsMat<3,3> Asch1 = {
		-149, -50, -154,
    	 537, 180,  546,
    	 -27,  -9,  -25,
	};	
	dispf(Asch1, "% 10.4f");
	ArcsMat<3,3> Phes1, Hhes1;
	Hessenberg(Asch1, Phes1, Hhes1);			// Hessenberg分解を計算 (引数渡し版)
	std::tie(Phes1, Hhes1) = Hessenberg(Asch1);	// Hessenberg分解を計算 (タプル返し版)
	dispf(Phes1, "% 8.4f");
	dispf(Hhes1, "% 10.4f");
	dispf(Phes1*Hhes1*~Phes1, "% 10.4f");		// 元に戻るかチェック
	constexpr auto PHhes2 = Hessenberg(Aqr1);	// コンパイル時にHessenberg分解を計算
	dispf(std::get<0>(PHhes2), "% 10.4f");		// コンパイル時に計算したユニタリ行列Pを表示
	dispf(std::get<1>(PHhes2), "% 10.4f");		// コンパイル時に計算したヘッセンベルグ行列Hを表示
	auto [Phes3, Hhes3] = Hessenberg(Acomp1);	// 複素数hessenberg分解を計算
	dispf(Phes3, "% 8.4f");
	dispf(Hhes3, "% 8.4f");
	dispf(Phes3*Hhes3*~Phes3, "% 8.4f");		// 元に戻るかチェック
	
	// Schur分解関連の関数
	printf("\n★★★★★★★ Schur分解関連の関数\n");
	dispf(Asch1, "% 10.4f");
	ArcsMat<3,3> Usch1, Tsch1;
	Schur(Asch1, Usch1, Tsch1);				// Schur分解を計算 (引数渡し版)
	std::tie(Usch1, Tsch1) = Schur(Asch1);	// Schur分解を計算 (戻り値返し版)
	dispf(Usch1, "% 10.4f");
	dispf(Tsch1, "% 10.4f");
	dispf(Usch1*Tsch1*~Usch1, "% 10.4f");	// 元に戻るかチェック
	dispf(~Usch1*Usch1, "% 10.4f");			// ユニタリ行列かチェック
	ArcsMat<3,3, std::complex<double>> Asch2 = {
		 4, -8,  1,
		 1, -5,  1,
		-9,  8, -6
	};
	dispf(Asch2, "% 10.4f");
	const auto [Usch2, Tsch2] = Schur(Asch2);	// 重根がある場合
	dispf(Usch2, "% 10.4f");
	dispf(Tsch2, "% 10.4f");
	dispf(Usch2*Tsch2*~Usch2, "% 10.4f");	// 元に戻るかチェック
	dispf(~Usch2*Usch2, "% 10.4f");			// ユニタリ行列かチェック
	
	/*
	ArcsMat<3,3> W = {
    -2.039207261322640e+00,     7.106155096368991e+00,    -8.161617539766114e+02,
    -5.733687193064796e-03,    -9.607928087369952e-01,     5.052140520240906e+01,
                       0.0,    -1.242272970283136e-09,     6.532248830737331e-08
	};
	const auto [Qw, Rw] = QR(W);
	dispf(Qw, "% 16.15e");
	*/

	/*
	dispf(Qsch1*Usch1*inv(Qsch1), "% 8.4f");	// 元に戻るかチェック
	constexpr auto QUsch1x = Schur(Asch1);		// コンパイル時にSchur分解を計算 (実行時と符号が異なる場合がある)
	constexpr auto Qsch1x = std::get<0>(QUsch1x);	// Qを取り出す
	constexpr auto Usch1x = std::get<1>(QUsch1x);	// Uを取り出す
	dispf(Qsch1x, "% 8.4f");
	dispf(Usch1x, "% 8.4f");
	dispf(Qsch1x*Usch1x*inv(Qsch1x), "% 8.4f");	// 元に戻るかチェック
	*/
	
	/*
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
	
	// 行列指数関数のテスト
	printf("\n★★★★★★★ 行列指数関数のテスト\n");
	Matrix<3,3> Y = expm(A, 6);
	PrintMat(A);
	PrintMatrix(Y,"% 16.14e");
	printf("det(Y)     = % 16.14e\n", det(Y));
	printf("exp(tr(A)) = % 16.14e\n\n", exp(tr(A)));	// 公式通りに一致！
	PrintMatrix(integral_expm(A,100e-6,10,6),"% 16.14e");
	
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

