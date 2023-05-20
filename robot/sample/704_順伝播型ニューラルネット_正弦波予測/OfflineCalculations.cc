//! @file OfflineCalculations.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2021/12/25
//! @author Yokokura, Yuki
//!
//! @par オフライン計算用のメインコード
//! - 「make offline」でコンパイルすると，いつものARCS制御用のコードは走らずに，
//!    このソースコードのみが走るようになる。
//! - ARCSライブラリはもちろんそのままいつも通り使用可能。
//! - 従って，オフラインで何か計算をしたいときに，このソースコードに記述すれば良い。
//!
// Copyright (C) 2011-2021 Yokokura, Yuki
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
#include "FrameGraphics.hh"
#include "CuiPlot.hh"
#include "SingleLayerPerceptron.hh"
#include "NeuralNetParamDef.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	printf("Predictive FFNN Test Code\n");
	
	// ここにオフライン計算のコードを記述
	
	// パラメータ設定＆初期設定
	constexpr size_t N = 1;		// [ch] 入力データチャネル数
	constexpr size_t K = 1;		// [ch] 出力/教師データチャネル数
	constexpr size_t D = 100;	// [samples] 時系列データの長さ
	constexpr size_t W = 12;	// [samples] 入力ウィンドウ幅
	constexpr size_t M = 10;	// [-]  ミニバッチ数
	constexpr size_t E = 1000;	// [-]	エポック数
	constexpr NnInitTypes INI = NnInitTypes::HE;			// 初期化方法
	constexpr NnDescentTypes OP = NnDescentTypes::MOMENTUM;	// 勾配降下方法
	constexpr ActvFunc AF1 = ActvFunc::ReLU;				// 活性化関数 レイヤ1
	constexpr ActvFunc AF2 = ActvFunc::ReLU;				// 活性化関数 レイヤ2
	constexpr ActvFunc AF3 = ActvFunc::IDENTITY;			// 活性化関数 レイヤ3
	constexpr NnDropout DO1 = NnDropout::DISABLE;			// ドロップアウト レイヤ1
	constexpr NnDropout DO2 = NnDropout::DISABLE;			// ドロップアウト レイヤ2
	constexpr NnDropout DO3 = NnDropout::DISABLE;			// ドロップアウト レイヤ3
	constexpr size_t P2 = 8;		// [-] パーセプトロンユニット数 レイヤ2
	FrameGraphics FG(1200, 1000);	// グラフプロット用フレームグラフィックス
	
	// 模擬観測データの生成
	Matrix<1,D> i = Matrix<1,D>::ramp();	// [-] サンプル時刻ベクトル
	Matrix<1,D> v;							// [-] 観測データベクトル
	for(size_t k = 1; k <= D; ++k){
		v[k] = sin(2.0*M_PI*2.0*(i[k] - 1)/(D - 1));	// [-] 正弦波観測データの模擬生成
	}
	CuiPlot PlotObsrv(FG, 0, 0, 1200, 250);
	PlotObsrv.SetRanges(0, i[D], -1.5, 1.5);
	PlotObsrv.SetGridDivision(10, 6);
	PlotObsrv.SetGridLabelFormat("%3.0f", "% 2.1f");
	PlotObsrv.SetAxisLabels("Time [samples]", "Observation Data [-]");
	PlotObsrv.DrawAxis();
	PlotObsrv.Plot(i, v, CuiPlotTypes::PLOT_STAIRS, FGcolors::CYAN);
	//PrintMatrix(v, "% 6.4f");
	
	// 入力データと教師データの整形(Many-to-Oneタイプ)
	Matrix<M,W> U;	// 入力データ
	Matrix<M,1>	R;	// 教師データ
	for(size_t k = 1; k <= M; ++k){
		U.SetVerticalVec(v.GetVerticalVec<W>(1, k), k, 1);
	}
	R = tp(v.GetVerticalVec<M>(1, W + 1));
	PrintMatrix(U, "% 6.4f");
	PrintMatrix(R, "% 6.4f");
	printf("Range of training samples = 1 to %zd [samples]\n", W + M);
	
	// 順伝播ニューラルネットの生成
	SingleLayerPerceptron<W, W,  M, AF1, INI, OP, DO1> L1;
	SingleLayerPerceptron<W, P2, M, AF2, INI, OP, DO2> L2;
	SingleLayerPerceptron<P2, K, M, AF3, INI, OP, DO3> L3;
	
	// ニューラルネットの初期設定
	L1.SetGainOfMomentumSGD(0.001, 0.9);
	L2.SetGainOfMomentumSGD(0.001, 0.9);
	L3.SetGainOfMomentumSGD(0.001, 0.9);
	L1.SetDropoutRate(0.95);
	L2.SetDropoutRate(0.95);
	L3.SetDropoutRate(0.95);
	L1.InitWeight(W);
	L2.InitWeight(W);
	L3.InitWeight(P2);
	
	// 誤差逆伝播用の変数
	Matrix<M,W>  Y1;
	Matrix<M,P2> Y2;
	Matrix<M,K>  Y;
	Matrix<M,P2> dY2;
	Matrix<M,W>  dY1;
	Matrix<M,W>  dU;
	
	// 誤差逆伝播法による訓練
	for(size_t e = 0; e <= E; ++e){
		// ドロップアウトの準備
		L1.CalcDropout();
		L2.CalcDropout();
		
		// 順伝播計算
		L1.CalcForwardForTraining(U, Y1);
		L2.CalcForwardForTraining(Y1, Y2);
		L3.CalcForwardForTraining(Y2, Y);
		
		// 逆伝播による誤差計算
		L3.CalcDeltaForOutputLayer(Y, R, dY2);
		L2.CalcDelta(dY2, dY1);
		L1.CalcDelta(dY1, dU);
		
		// 重みとバイアスの更新
		L1.UpdateWeight(U);
		L2.UpdateWeight(Y1);
		L3.UpdateWeight(Y2);
		
		printf("%5zu : %f\n", e, L3.GetLoss(Y, R));
	}
	
	// プロット出力
	FG.SavePngImageFile("FFNNPlot.png");
	
	return EXIT_SUCCESS;	// 正常終了
}

