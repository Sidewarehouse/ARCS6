//! @file OfflineCalculations.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2022/07/22
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

// 追加のARCSライブラリをここに記述
#include "Matrix.hh"
#include "CsvManipulator.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	// 行列宣言とセットのテスト
	printf("◆ 行列宣言とセットのテスト\n");
	Matrix<3,3> A = {	// 宣言と同時に値をセットする場合
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	};
	Matrix<3,3> B;		// 宣言したあとに、
	B.Set(				// 値をセットする場合
		5,  1, -3,
		3, -1,  7,
		4,  2, -5
	);
	auto C = A;			// 宣言と同時に既にある行列を代入する場合
	
	// 行列の表示のテスト
	printf("\n◆ 行列の表示のテスト\n");
	PrintMat(A);			// 行列の表示
	PrintMat(B);
	PrintMatrix(C,"% 6.4f");// 表示の書式指定をする場合
	PrintMatSize(A);		// 行列のサイズの表示
	
	// 定数行列，コンパイル時定数行列と初期化のテスト
	printf("\n◆ 定数行列，コンパイル時定数行列と初期化のテスト\n");
	const Matrix<2,2> Alpha = {
		1, 2,
		3, 4
	};
	constexpr Matrix<2,2> Beta = {
		5, 6,
		7, 8
	};
	PrintMat(Alpha);
	PrintMat(Beta);
	
	// 縦ベクトルの[]オペレータによる要素アクセスのテスト
	printf("\n◆ 縦ベクトルの[]オペレータによる要素アクセスのテスト\n");
	Matrix<1,5> alpha = {
		1,
		2,
		3,
		4,
		5
	};
	PrintMat(alpha);
	printf("alpha[3] = % g\n", alpha[3]);
	alpha[3] = 6;
	printf("alpha[3] = % g\n\n", alpha[3]);
	
	// 行列の()オペレータによる要素アクセスのテスト
	printf("\n◆ 行列の()オペレータによる要素アクセスのテスト\n");
	printf("B(1,1) = % g\n", B(1,1));
	printf("B(1,2) = % g\n", B(1,2));
	printf("B(1,3) = % g\n", B(1,3));
	printf("B(2,1) = % g\n", B(2,1));
	printf("B(2,2) = % g\n", B(2,2));
	printf("B(2,3) = % g\n", B(2,3));
	printf("B(3,1) = % g\n", B(3,1));
	printf("B(3,2) = % g\n", B(3,2));
	printf("B(3,3) = % g\n", B(3,3));
	B(2,3) = 0;
	PrintMat(B);
	B(2,3) = 7;
	PrintMat(B);
	
	// 演算子のテスト
	printf("\n◆ 演算子のテスト\n");
	C = +A;
	printf("+A = "); PrintMat(C);
	C = -A;
	printf("-A = "); PrintMat(C);
	C = A + B;
	printf("A + B = "); PrintMat(C);
	C = A + 3;
	printf("A + 3 = "); PrintMat(C);
	C = A - 3;
	printf("A - 3 = "); PrintMat(C);
	C = A - B;
	printf("A - B = "); PrintMat(C);
	C = A*B;
	printf("A*B = "); PrintMat(C);
	C = A*2;
	printf("A*2 = "); PrintMat(C);
	C += A;
	printf("+=A = "); PrintMat(C);
	C += 2;
	printf("+=2 = "); PrintMat(C);
	C -= A;
	printf("-=A = "); PrintMat(C);
	C -= 1;
	printf("-=1 = "); PrintMat(C);
	C = A^2;
	printf("A^2 = "); PrintMat(C);
	C = 3 + A;
	printf("3 + A = "); PrintMat(C);
	C = 3 - A;
	printf("3 - A = "); PrintMat(C);
	C = 2*A;
	printf("2*A = "); PrintMat(C);
	
	// アダマール乗除演算子のテスト
	printf("\n◆ アダマール乗除演算子のテスト\n");
	auto Gamma = Alpha & Beta;	// アダマール積
	PrintMat(Gamma);
	Gamma = Alpha % Beta;		// アダマール除算
	PrintMat(Gamma);
	
	// 行列の生配列のテスト
	printf("\n◆ 行列の生配列のテスト\n");
	Matrix<3,3> J[3];
	for(size_t i = 0; i < 3; ++i){
		J[i].Set(
			i+1,  i+1,  i+1,
			  2,    3,   -2,
			  3,   -1,  i+1
		);
		PrintMat(J[i]);
	}
	for(size_t i = 0; i < 3; ++i) PrintMatrix(inv(A*J[i]),"% 6.4f");	// 演算例
	
	// 行列のstd::array配列のテスト
	printf("\n◆ 行列のstd::array配列のテスト\n");
	std::array<Matrix<3,3>, 3> K;
	for(size_t i = 0; i < 3; ++i){
		K.at(i).Set(
			i+1,  i+1,  i+1,
			  2,    3,   -2,
			  3,   -1,  i+1
		);
		PrintMat(K.at(i));
	}
	for(size_t i = 0; i < 3; ++i) PrintMatrix(inv(A*K.at(i)),"% 6.4f");	// 演算例
	
	// 零行列，1行列，単位行列のテスト
	printf("\n◆ 零行列，1行列，単位行列のテスト\n");
	auto O = Matrix<3,3>::zeros();
	auto l = Matrix<3,3>::ones();
	constexpr auto I = Matrix<3,3>::eye();
	PrintMat(O);
	PrintMat(l);
	PrintMat(I);
	
	// 行列要素に設定＆から取得する系の関数のテスト
	printf("\n◆ 行列要素に設定＆から取得する系の関数のテスト\n");
	Matrix<2,3> D;
	D.Set(
		-1,2,
		-3,4,
		-5,6
	);
	auto E = D;
	Matrix<2,1> w;
	w.Set(
		11,22
	);
	Matrix<1,3> v;
	v.Set(
		11,
		22,
		33
	);
	PrintMat(D);
	double d11, d12, d21, d22, d31, d32;
	D.Get(d11, d12, d21, d22, d31, d32);
	printf("d11 = %f, d12 = %f, d21 = %f, d22 = %f, d31 = %f, d32 = %f\n\n", d11, d12, d21, d22, d31, d32);
	std::array<double, 10> v10array = {100,200,300,400,500,600,700,800,900,1000};
	Matrix<1,10> v10;
	v10.LoadArray(v10array);
	v10.LoadShortVector(v);
	PrintMat(v10);
	D.SetElem(3,1,3.14159);
	PrintMat(D);
	printf("D.GetElem = %f\n\n", D.GetElem(3,1));
	D.SetElement(1,3,2.71828);
	PrintMat(D);
	printf("D.GetElement = %f\n\n", D.GetElement(1,3));
	PrintMat(v10.GetVerticalVec<3>(1,3));
	v10.SetVerticalVec(v, 1, 4);
	PrintMat(v10);
	PrintMat(B.GetHorizontalVec<2>(2,3));
	B.SetHorizontalVec(tp(v), 1, 1);
	PrintMat(B);
	E.FillAll(3.939);
	PrintMat(E);
	E.FillAllZero();
	PrintMat(E);
	
	// 行列演算補助系の関数のテスト
	printf("\n◆ 行列演算補助系の関数のテスト\n");
	PrintMat(tp(D));
	printf("tr(A) = %f\n", tr(A));
	printf("prod(A) = %f\n", prod(A));
	PrintMat(diag(A));
	PrintMat(sumrow(D));
	PrintMat(sumcolumn(D));
	printf("max(v) = %f\n", max(v));
	printf("absmax(v) = %f\n", absmax(v));
	printf("maxidx(v) = %ld\n", maxidx(v));
	printf("absmaxidx(v) = %ld\n", absmaxidx(v));
	printf("nonzeroele(I) = %ld\n", nonzeroele(I));
	printf("rank(I) = %ld\n", rank(I));
	
	// 要素ごとの数学関数のテスト
	printf("\n◆ 要素ごとの数学関数のテスト\n");
	PrintMat( expe( Matrix<3,3>::eye())         );
	PrintMat( loge( Matrix<3,3>::eye()*2.71828) );
	PrintMat( abse(-Matrix<3,3>::eye())         );
	PrintMat( sqrte(Matrix<3,3>::eye()*2.0)     );
	
	// 行・列操作系関数のテスト
	printf("\n◆ 行・列操作系の関数のテスト\n");
	PrintMat(getrow(D, 2));
	setrow(D, w, 2);
	PrintMat(D);
	swaprow(D, 1, 2);
	PrintMat(D);
	fillrow(D, 3.14, 1,1,2);
	PrintMat(D);
	PrintMat(getcolumn(D, 1));
	setcolumn(D, v, 1);
	PrintMat(D);
	swapcolumn(D, 1, 2);
	PrintMat(D);
	fillcolumn(D, 3.14, 1, 1, 2);
	PrintMat(D);
	Matrix<6,6> D2 = {
		 1,  2,  3,  4,  5,  6,
		 7,  8,  9, 10, 11, 12,
		13, 14, 15, 16, 17, 18,
		19, 20, 21, 22, 23, 24,
		25, 26, 27, 28, 29, 30,
		31, 32, 33, 34, 35, 36
	};
	setvvector(D2, getcolumn(D, 1), 2, 4);
	PrintMatrix(D2, "% 6.2f");
	Matrix<1,3> vD2;
	getvvector(D2, 2, 4, vD2);
	PrintMat(vD2);
	Matrix<3,1> hD2;
	gethvector(D2, 5, 3, hD2);
	PrintMat(hD2);
	Matrix<3,3> sD2;
	getsubmatrix(D2, 2, 4, sD2);
	PrintMat(sD2);
	setsubmatrix(D2, 2, 3, A);
	PrintMat(D2);
	
	// 並び替え系関数のテスト
	printf("\n◆ 並び替え系の関数のテスト\n");
	Matrix<1,3,int> vi = {
		3,
		1,
		2
	};
	PrintMat(orderrow(D, vi));
	PrintMat(reorderrow(D, vi));
	PrintMat(shiftup(A));
	PrintMat(shiftdown(A));
	PrintMat(shiftright(A));
	PrintMat(shiftleft(A));
	PrintMat(shiftup(A, 2));
	PrintMat(shiftdown(A, 2));
	PrintMat(shiftright(A, 2));
	PrintMat(shiftleft(A, 2));
	PrintMat(gettriup(A));
	PrintMat(gettriup(A, 1));
	
	// ノルム演算系のテスト
	printf("\n◆ ノルム演算系のテスト\n");
	printf("infnorm(A) = %f\n", infnorm(A));
	printf("euclidnorm(v) = %f\n", euclidnorm(v));
	
	// LU分解のテスト
	printf("\n◆ LU分解のテスト\n");
	Matrix<1,3,int> vv;			// 並べ替え記憶列ベクトル
	Matrix<3,3> L, U;
	LU(A,L,U,vv);				// LU分解
	PrintMat(L);
	PrintMat(U);
	PrintMat(L*U);					// 元の行列に戻るか確認。私達、入れ替わってる！
	PrintMat(reorderrow(L*U, vv));	// なので行の順番を元に戻す。入力行列と一緒！
	
	// Cholesky分解(LDL^T版)のテスト
	printf("\n◆ Cholesky分解(LDL^T版)のテスト\n");
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
	printf("\n◆ Cholesky分解(LL^T版)のテスト\n");
	Cholesky(Ach, Lch);
	PrintMat(Ach);
	PrintMat(Lch);
	PrintMat(Lch*tp(Lch));
	
	// QR分解のテスト1
	printf("\n◆ QR分解(実数版)のテスト1\n");
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
	printf("\n◆ QR分解(実数版)のテスト2\n");
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
	printf("\n◆ QR分解(実数版)のテスト3\n");
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
	printf("\n◆ QR分解(実数版)のテスト4\n");
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
	printf("\n◆ SVD特異値分解のテスト1(縦長行列の場合)\n");
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
	printf("\n◆ SVD特異値分解のテスト2(横長行列の場合)\n");
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
	printf("\n◆ SVD特異値分解のテスト3(ランク落ちの場合)\n");
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
	printf("\n◆ SVD特異値分解のテスト4(符号修正が必要な場合)\n");
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
	printf("\n◆ Schur分解のテスト1(実数固有値の場合)\n");
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
	printf("\n◆ Schur分解のテスト2(複素数固有値の場合)\n");
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
	printf("\n◆ 連立方程式の球解テスト\n");
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
	printf("\n◆ 上三角行列の連立方程式の球解テスト\n");
	Matrix<3,3> Auptri = {
		1, 3, 6,
		0, 2, 7,
		0, 0,-4
	};
	solve_upper_tri(Auptri, b, xslv);
	PrintMat(Auptri);
	PrintMat(xslv);
	
	// 行列式と逆行列のテスト
	printf("\n◆ 行列式と逆行列のテスト\n");
	printf("det(A) = %f\n", det(A));
	PrintMatrix(inv(A), "% 16.14e");
	PrintMatrix(inv_with_check(A), "% 16.14e");
	
	// 左上小行列の逆行列のテスト
	printf("\n◆ 左上小行列の逆行列のテスト\n");
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
	printf("\n◆ 上三角行列の逆行列のテスト\n");
	Matrix<3,3> Auptri_inv;
	inv_upper_tri(Auptri, Auptri_inv);
	PrintMat(Auptri_inv);
	
	// 上三角行列で左上小行列の逆行列のテスト
	printf("\n◆ 上三角行列で左上小行列の逆行列のテスト\n");
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
	printf("\n◆ 疑似逆行列のテスト\n");
	PrintMatrix(lpinv(D), "% 16.14e");
	PrintMatrix(rpinv(tp(D)), "% 16.14e");
	Matrix<1,2> Dpinv = {
		1,
		2
	};
	PrintMatrix(lpinv(Dpinv), "% 16.14e");
	PrintMatrix(rpinv(tp(Dpinv)), "% 16.14e");
	
	// 左上小行列の疑似逆行列のテスト
	printf("\n◆ 左上小行列の疑似逆行列のテスト\n");
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
	printf("\n◆ 行列指数関数のテスト\n");
	Matrix<3,3> Y = expm(A, 6);
	PrintMat(A);
	PrintMatrix(Y,"% 16.14e");
	printf("det(Y)     = % 16.14e\n", det(Y));
	printf("exp(tr(A)) = % 16.14e\n\n", exp(tr(A)));	// 公式通りに一致！
	PrintMatrix(integral_expm(A,100e-6,10,6),"% 16.14e");
	
	// 特に意味のないリッカチ方程式のテスト
	printf("\n◆ 特に意味のないリッカチ方程式のテスト\n");
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
	printf("◆ float型のテスト\n");
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
	printf("◆ int型のテスト\n");
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
	printf("◆ long型のテスト\n");
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
	printf("◆ 複素数型のテスト\n");
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
	printf("◆ 負の平方根のテスト\n");
	Matrix<3,3,std::complex<double>> Acomp = {
		std::complex(1.0,0.0), std::complex( 1.0,0.0), std::complex( 1.0,0.0),
		std::complex(2.0,0.0), std::complex( 3.0,0.0), std::complex(-2.0,0.0),
		std::complex(3.0,0.0), std::complex(-1.0,0.0), std::complex( 1.0,0.0),
	};
	PrintMat(Acomp);
	PrintMatrix(sqrte(A), "% 6.3f");	// double型だとnanになってしまうが，
	PrintMatrix(sqrte(Acomp), "% 6.3f");// std::complexだとちゃんと計算できる
	
	// エルミート転置のテスト
	printf("◆ エルミート転置のテスト\n");
	Matrix<2,3,std::complex<double>> Acomp2 = {
		std::complex(1.0, 2.0), std::complex(  3.0, -4.0),
		std::complex(5.0,-6.0), std::complex(  7.0,  8.0),
		std::complex(9.0,10.0), std::complex(-11.0,-12.0)
	};
	PrintMat(Acomp2);
	PrintMat(Htp(Acomp2));
	
	// 複素数LU分解のテスト
	printf("◆ 複素数LU分解のテスト\n");
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
	printf("◆ 複素数逆行列のテスト\n");
	PrintMat(inv(Acomp3));
	PrintMatrix(inv(Acomp3)*Acomp3, "% 7.3f");
	
	// 複素数QR分解のテスト
	printf("◆ 複素数QR分解のテスト\n");
	Matrix<3,3,std::complex<double>> Qcqr, Rcqr;
	QR(Acomp3, Qcqr, Rcqr);
	PrintMatrix(Qcqr, "% 7.3f");
	PrintMatrix(Rcqr, "% 7.3f");
	PrintMatrix(Qcqr*Htp(Qcqr), "% 7.3f");	// Qが直交行列かチェック
	PrintMat(Qcqr*Rcqr);					// 元に戻るかチェック
	
	// 固有値計算のテスト1(実数固有値の場合)
	printf("◆ 固有値計算のテスト1(実数固有値の場合)\n");
	Matrix<3,3> Aeig = {
		-3, -4,  2,
		-7,  1, -5,
		 6, -7,  3
	};
	PrintMat(Aeig);
	PrintMat(eigen(Aeig));
	PrintMat(eigenvec(Aeig));
	
	// 固有値計算のテスト2(複素数固有値の場合)
	printf("◆ 固有値計算のテスト2(複素数固有値の場合)\n");
	Aeig.Set(
		10,  -8,   5,
		-8,   9,   6,
		-1, -10,   7
	);
	PrintMat(Aeig);
	PrintMat(eigen(Aeig));
	PrintMat(eigenvec(Aeig));
	
	// 2乗のコンパイル時定数演算のテスト
	printf("\n◆ 2乗のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cx = {
		1,  1,  1,
		2,  3, -2,
		3, -1,  1
	};
	constexpr Matrix<3,3> Cxsq = Cx^2;			// 2乗のコンパイル時計算
	PrintMat(Cxsq);
	
	// 逆行列のコンパイル時定数演算のテスト
	printf("\n◆ 逆行列のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cxinv = inv(Cx);		// 逆行列のコンパイル時計算
	PrintMat(Cxinv);
	
	// 状態遷移行列のコンパイル時定数演算のテスト
	printf("\n◆ 状態遷移行列のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cxexp = expm(Cx, 6);	// 状態遷移行列のコンパイル時計算
	PrintMatrix(Cxexp, "% 16.14e");
	printf("det(Y)     = % 16.14e\n", det(Cxexp));
	printf("exp(tr(A)) = % 16.14e\n\n", exp(tr(Cx)));	// 公式通りに一致！
	
	// 状態遷移行列の定積分のコンパイル時定数演算のテスト
	printf("\n◆ 状態遷移行列の定積分のコンパイル時定数演算のテスト\n");
	constexpr Matrix<3,3> Cxexpint = integral_expm(Cx, 100e-6, 10, 6);	// 状態遷移行列の定積分のコンパイル時計算
	PrintMatrix(Cxexpint, "% 16.14e");
	
	// SVD特異値分解のコンパイル時定数演算のテスト
	printf("\n◆ SVD特異値分解のコンパイル時定数演算のテスト\n");
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
	printf("\n◆ 行列のランクのコンパイル時定数演算のテスト\n");
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
	printf("\n◆ クロネッカー積のテスト\n");
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
	printf("\n◆ vec作用素のテスト\n");
	PrintMat(vec(A));
	PrintMat(vec(Axsvd));
	Matrix<1,9>::vecinv(vec(A), A);			// 引数で返す版
	A = Matrix<1,9>::vecinv<3,3>(vec(A));	// 戻り値で返す版
	PrintMat(A);
	auto Avec = Matrix<1,8>::vecinv<2,4>(vec(Axsvd));	// 戻り値で返す版
	PrintMat(Avec);
	
	return EXIT_SUCCESS;	// 正常終了
}

