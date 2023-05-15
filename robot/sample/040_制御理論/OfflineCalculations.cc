//! @file OfflineCalculations.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2022/07/27
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
#include "ArcsControl.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	
	// プラント状態空間モデル(適当)
	printf("◆ プラント状態空間モデル(適当)\n");
	constexpr Matrix<3,3> A = {
		-1, -2, -9,
    	 0, -1, -7,
    	-2,  6, -1
	};
	constexpr Matrix<1,3> b = {1, 0, 0};
	constexpr Matrix<3,1> c = {1, 0, 0};
	PrintMat(A);
	PrintMat(b);
	PrintMat(c);
	
	
	printf("◆ 連続リアプノフ方程式の解(MATLABでいうlyap)のテスト\n");
	constexpr auto Q = Matrix<3,3>::eye();
	Matrix<3,3> X;
	ArcsControl::Lyapunov(A, Q, X);		// A*X + X*A^T + Q = 0 の解X(引数で返す版)
	X = ArcsControl::Lyapunov(A, Q);	// A*X + X*A^T + Q = 0 の解X(戻り値で返す版)
	PrintMat(X);
	PrintMat(A*X + X*tp(A) + Q);					// 零行列になればOK
	
	constexpr auto Xx = ArcsControl::Lyapunov(A, Q);// コンパイル時定数版
	PrintMat(Xx);
	
	
	printf("◆ 可制御/可観測グラミアン(MATLABでいうgram)のテスト\n");
	Matrix<3,3> Wc, Wo;
	ArcsControl::GramianCtrl(A, b, Wc);		// 可制御グラミアン(引数で返す版)
	Wc = ArcsControl::GramianCtrl(A, b);	// 可制御グラミアン(戻り値で返す版)
	PrintMat(Wc);
	
	ArcsControl::GramianObsrv(A, c, Wo);	// 可観測グラミアン(引数で返す版)
	Wo = ArcsControl::GramianObsrv(A, c);	// 可観測グラミアン(戻り値で返す版)
	PrintMat(Wo);
	
	constexpr auto Wcx = ArcsControl::GramianCtrl(A, b);	// コンパイル時定数版
	PrintMat(Wcx);
	constexpr auto Wox = ArcsControl::GramianObsrv(A, c);	// コンパイル時定数版
	PrintMat(Wox);
	
	
	printf("◆ 平衡実現(入出力平衡化、MATLABでいうbalreal)のテスト\n");
	Matrix<3,3> Ah;
	Matrix<1,3> bh;
	Matrix<3,1> ch;
	ArcsControl::BalanceReal(A, b, c, Ah, bh, ch);				// 引数で返す版
	std::tie(Ah, bh, ch) = ArcsControl::BalanceReal(A, b, c);	// タプルで返す版
	PrintMat(Ah);
	PrintMat(bh);
	PrintMat(ch);
	auto O = ArcsControl::GramianCtrl(Ah, bh) - ArcsControl::GramianObsrv(Ah, ch);
	PrintMatrix(O, "% 6.4f");	// 平衡化後は可制御と可観測グラミアンが一緒になる
	
	constexpr auto Abchx = ArcsControl::BalanceReal(A, b, c);	// コンパイル時定数版
	PrintMat(std::get<0>(Abchx));	// constexprの初期化に構造化束縛は使えない。
	PrintMat(std::get<1>(Abchx));	// そこで、std::get<n>で順に取り出す。
	PrintMat(std::get<2>(Abchx));
	
	
	printf("◆ 離散化(MATLABでいうc2d)のテスト\n");
	// 2慣性共振系のパラメータ例
	constexpr double Ts = 100e-6; // [s]      サンプリング時間
	constexpr double Jm = 1e-4;   // [kgm^2]  モータ慣性
	constexpr double Jl = 1;      // [kgm^2]  負荷側慣性
	constexpr double Dm = 0.1;    // [Nms/rad]モータ粘性
	constexpr double Dl = 0.3;    // [Nms/rad]負荷側粘性
	constexpr double Ks = 500;    // [Nm/rad] 2慣性間のばね定数
	constexpr double Rg = 50;     //          減速比
	constexpr double Kt = 0.04;   // [Nm/A]   トルク定数
	constexpr Matrix<3,3> Ac = {
		-Dm/Jm,      0, -Ks/(Rg*Jm),
			 0, -Dl/Jl,       Ks/Jl,
		1.0/Rg,     -1,           0
	};
	constexpr Matrix<2,3> Bc = {
		Kt/Jm,       0,
			0, -1.0/Jl,
			0,       0
	};
	constexpr Matrix<3,1> C = {1, 0, 0};
	
	auto [Ad, Bd] = ArcsControl::Discretize(Ac, Bc, Ts);// タプルで返す版
	ArcsControl::Discretize(Ac, Bc, Ad, Bd, Ts);		// 引数で返す版
	PrintMat(Ad);
	PrintMat(Bd);
	
	std::tie(Ad, Bd) = ArcsControl::Discretize(Ac, Bc, Ts, 1, 1);	// パデ近似の次数と積分分割数を指定する版
	PrintMat(Ad);
	PrintMat(Bd);
	
	constexpr auto AdBd = ArcsControl::Discretize(Ac, Bc, Ts, 1, 1);// コンパイル時定数版
	PrintMat(std::get<0>(AdBd));	// constexprの初期化に構造化束縛は使えない。
	PrintMat(std::get<1>(AdBd));	// そこで、std::get<n>で順に取り出す。
	
	auto [Ach, Bch, Ch] = ArcsControl::BalanceReal(Ac, Bc, C);	// 平衡化してから、
	std::tie(Ad, Bd) = ArcsControl::Discretize(Ach, Bch, Ts);	// 離散化する場合
	PrintMat(Ad);
	PrintMat(Bd);
	
	return EXIT_SUCCESS;	// 正常終了
}

