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
#include "ArcsControl.hh"

using namespace ARCS;
using namespace ArcsMatrix;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	
	// ここにオフライン計算のコードを記述
	// プラント状態空間モデル(適当)
	printf("◆ プラント状態空間モデル(適当)\n");
	constexpr ArcsMat<3,3> Ap = {
		-1, -2, -9,
    	 0, -1, -7,
    	-2,  6, -1
	};
	constexpr ArcsMat<3,1> bp = {
		1,
		0,
		0};
	constexpr ArcsMat<1,3> cp = {1, 0, 0};
	disp(Ap);
	disp(bp);
	disp(cp);

	printf("◆ 連続リアプノフ方程式の解(MATLABでいうlyap)\n");
	constexpr auto Qp = ArcsMat<3,3>::eye();
	ArcsMat<3,3> Xp;
	ArcsControl::Lyapunov(Ap, Qp, Xp);		// AX + XA' + Q = 0 の解X (引数渡し版)
	Xp = ArcsControl::Lyapunov(Ap, Qp);		// AX + XA' + Q = 0 の解X (戻り値返し版)
	dispf(Xp, "% 8.4f");
	constexpr auto Xpx = ArcsControl::Lyapunov(Ap, Qp);		// コンパイル時に連続リアプノフ方程式を計算
	dispf(Xpx, "% 8.4f");
	printf("|| AX + XA' + Q || = %e\n\n", norm(Ap*Xp + Xp*~Ap + Qp));	// ←ほぼゼロならオーケー

	printf("◆ 離散リアプノフ方程式の解(MATLABでいうdlyap)\n");
	ArcsControl::DiscLyapunov(Ap, Qp, Xp);	// AXA' - X + Q = 0 の解X (引数渡し版)
	Xp = ArcsControl::DiscLyapunov(Ap, Qp);	// AXA' - X + Q = 0 の解X (戻り値返し版)
	dispf(Xp, "% 8.4f");
	constexpr auto Xpx2 = ArcsControl::DiscLyapunov(Ap, Qp);// コンパイル時に離散リアプノフ方程式を解く
	dispf(Xpx2, "% 8.4f");
	printf("|| AXA' - X + Q || = %e\n\n", norm(Ap*Xp*~Ap - Xp + Qp));	// ←ほぼゼロならオーケー

	printf("◆ 原点零点がある場合のリアプノフ方程式の解\n");
	ArcsMat<2,2> Az = {
		-3.4182, 0.0,
		 1.0,    0.0
	};
	ArcsMat<2,2> Bz = {
		4.4545e3, -9.0909e3,
		0.0,       0.0
	};
	dispf(ArcsControl::Lyapunov(Az, Bz*~Bz), "% 8.4f");	// リアプノフ方程式の一意な解が無い。NaNが出現。

	printf("◆ 可制御/可観測グラミアン(MATLABでいうgram)\n");
	ArcsMat<3,3> Wpc, Wpo;
	ArcsControl::GramianCtrl(Ap, bp, Wpc);	// 可制御グラミアン (引数渡し版)
	Wpc = ArcsControl::GramianCtrl(Ap, bp);	// 可制御グラミアン (戻り値返し版)
	dispf(Wpc, "% 8.4f");
	constexpr auto Wpcx = ArcsControl::GramianCtrl(Ap, bp);	// コンパイル時に可制御グラミアンを計算
	dispf(Wpcx, "% 8.4f");
	ArcsControl::GramianObsrv(Ap, cp, Wpo);	// 可観測グラミアン (引数渡し版)
	Wpo = ArcsControl::GramianObsrv(Ap, cp);// 可観測グラミアン (戻り値返し版)
	dispf(Wpo, "% 8.4f");
	constexpr auto Wpox = ArcsControl::GramianObsrv(Ap, cp);// コンパイル時に可観測グラミアンを計算
	dispf(Wpox, "% 8.4f");

	printf("◆ 平衡実現(入出力平衡化、MATLABでいうbalreal)\n");
	ArcsMat<3,3> Aph;
	ArcsMat<3,1> bph;
	ArcsMat<1,3> cph;
	ArcsControl::BalanceReal(Ap, bp, cp, Aph, bph, cph);			// 平衡実現 (引数渡し版)
	std::tie(Aph, bph, cph) = ArcsControl::BalanceReal(Ap, bp, cp);	// 平衡実現 (タプル返し版)
	dispf(Aph, "% 8.4f");
	dispf(bph, "% 8.4f");
	dispf(cph, "% 8.4f");
	constexpr auto Abcph = ArcsControl::BalanceReal(Ap, bp, cp);	// コンパイル時に平衡実現
	dispf(std::get<0>(Abcph), "% 8.4f");
	dispf(std::get<1>(Abcph), "% 8.4f");
	dispf(std::get<2>(Abcph), "% 8.4f");
	const double gram1 = norm(ArcsControl::GramianCtrl(Ap, bp) - ArcsControl::GramianObsrv(Ap, cp));
	const double gram2 = norm(ArcsControl::GramianCtrl(Aph, bph) - ArcsControl::GramianObsrv(Aph, cph));
	printf("|| Wc - Wo || = %e\n",   gram1);	// 平衡化前は可制御グラミアンと可観測グラミアンは違うが、
	printf("|| Wc - Wo || = %e\n\n", gram2);	// 平衡化後は可制御グラミアンと可観測グラミアンが一緒になる
	
	printf("◆ 離散化(MATLABでいうc2d)\n");
	// 2慣性共振系のパラメータ例
	constexpr double Ts = 100e-6;	// [s]      サンプリング時間
	constexpr double Jl = 1;		// [kgm^2]  負荷側慣性
	constexpr double Dl = 0.3;		// [Nms/rad]負荷側粘性
	constexpr double Ds = 1;		// [Nms/rad]ねじれ粘性
	constexpr double Ks = 500;		// [Nm/rad] 2慣性間のばね定数
	constexpr double Jm = 1e-4;		// [kgm^2]  モータ慣性
	constexpr double Dm = 0.1;		// [Nms/rad]モータ粘性
	constexpr double Rg = 100;		//          減速比
	constexpr double Kt = 0.04;		// [Nm/A]   トルク定数
	// 状態ベクトルの定義 xp = [ωl, θs, ωm]^T
	constexpr ArcsMat<3,3> Atc = {
		-(Dl + Ds)/Jl,	Ks/Jl,			Ds/(Jl*Rg)			  ,
		-1,				0,				1.0/Rg				  ,
		Ds/(Jm*Rg),		-Ks/(Jm*Rg),	-(Ds/(Rg*Rg) + Dm)/Jm
	};
	constexpr ArcsMat<3,2> Btc = {
		0		, -1.0/Jl,
		0		, 0	     ,
		Kt/Jm	, 0
	};
	constexpr ArcsMat<3,1> Btc2 = {0, 0, Kt/Jm};
	constexpr auto Ct = ArcsMat<3,3>::eye();
	constexpr ArcsMat<1,3> Ct2 = {0, 0, 1};
	ArcsMat<3,3> Atd;
	ArcsMat<3,2> Btd;
	ArcsControl::Discretize(Atc, Btc, Atd, Btd, Ts);			// 離散化 (引数渡し版)
	std::tie(Atd, Btd) = ArcsControl::Discretize(Atc, Btc, Ts);	// 離散化 (タプル返し版)
	dispf(Atd, "% 8.4f");
	dispf(Btd, "% 8.4f");
	constexpr auto ABtd = ArcsControl::Discretize<13,30>(Atc, Btc, Ts);	// コンパイル時に離散化
	dispf(std::get<0>(ABtd), "% 8.4f");	// ↑ パデ近似の次数を13, 積分分割数を30に設定
	dispf(std::get<1>(ABtd), "% 8.4f");	// 積分分割数を増やしすぎると、‘constexpr’ evaluation operation count exceeds limit of 33554432 エラーが出る

	printf("◆ システムの判定\n");
	const bool StblAtc = ArcsControl::IsStable(Atc);			// 安定性判別 (MATLABでいうisstable)
	printf("IsStable(Atc) = %s\n", (StblAtc ? "true" : "false"));
	const bool StblAp = ArcsControl::IsStable(Ap);				// 安定性判別
	printf("IsStable(Ap)  = %s\n", (StblAp ? "true" : "false"));
	const bool StblAz = ArcsControl::IsStable(Az);				// 安定性判別
	printf("IsStable(Az)  = %s\n", (StblAz ? "true" : "false"));
	ArcsMat<2,2> Ao = {
		1,  1,
		4, -2
	};
	ArcsMat<2,2> Bo = {
		1, -1,
		1, -1
	};
	ArcsMat<2,2> Co = {
		-1,  1,
		 1, -1
	};
	const bool ObsvAo   = ArcsControl::IsObsv(Ao, Co);
	printf("IsObsv(Ao, Co)   = %s\n", (ObsvAo ? "true" : "false"));		// 可観測性判定 (MATLABでいう length(A) == rank(obsv(A,C))
	const bool ObsvAtc  = ArcsControl::IsObsv(Atc, Ct);
	printf("IsObsv(Atc, Ct)  = %s\n", (ObsvAtc ? "true" : "false"));	// 可観測性判定
	const bool ObsvAtc2 = ArcsControl::IsObsv(Atc, Ct2);
	printf("IsObsv(Atc, Ct2) = %s\n", (ObsvAtc2 ? "true" : "false"));	// 可観測性判定
	const bool CtrbAo   = ArcsControl::IsCtrb(Ao, Bo);
	printf("IsCtrb(Ao, Bo)    = %s\n", (CtrbAo ? "true" : "false"));	// 可制御性判定 (MATLABでいう length(A) == rank(ctrb(A,B))
	const bool CtrbAtc  = ArcsControl::IsCtrb(Atc, Btc);
	printf("IsCtrb(Atc, Btc)  = %s\n", (CtrbAtc ? "true" : "false"));	// 可制御性判定
	const bool CtrbAtc2 = ArcsControl::IsCtrb(Atc, Btc2);
	printf("IsCtrb(Atc, Btc2) = %s\n", (CtrbAtc2 ? "true" : "false"));	// 可制御性判定

	ArcsControl::DiscStateSpace<3, 1, 1> Sys1;

	return EXIT_SUCCESS;	// 正常終了
}

