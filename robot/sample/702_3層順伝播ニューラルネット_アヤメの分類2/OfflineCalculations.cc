//! @file OfflineCalculations.cc
//! @brief ARCS6 オフライン計算用メインコード
//! @date 2020/05/17
//! @author Yokokura, Yuki
//!
//! @par オフライン計算用のメインコード
//! - 「make offline」でコンパイルすると，いつものARCS制御用のコードは走らずに，
//!    このソースコードのみが走るようになる。
//! - ARCSライブラリはもちろんそのままいつも通り使用可能。
//! - 従って，オフラインで何か計算をしたいときに，このソースコードに記述すれば良い。
//!
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.
//
// 3層ニューラルネットワークのテスト用コード
//

// 基本のインクルードファイル
#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <array>

// 追加のARCSライブラリをここに記述
#include "Matrix.hh"
#include "CsvManipulator.hh"
#include "IrisClassDatasets.hh"
#include "FeedforwardNeuralNet3.hh"

using namespace ARCS;

//! @brief エントリポイント
//! @return 終了ステータス
int main(void){
	printf("ARCS OFFLINE CALCULATION MODE\n");
	printf("3-Layers Perceptron TEST CODE\n");
	
	// パラメータ設定
	constexpr size_t InSize  = 4;		// 入力層のサイズ
	constexpr size_t HidSize = 5;		// 隠れ層のサイズ
	constexpr size_t OutSize = 3;		// 出力層のサイズ 
	constexpr size_t MinbatSize = 10;	// ミニバッチサイズ
	constexpr size_t MinbatNum  = 12;	// ミニバッチ数
	constexpr size_t Epoch = 10000;		// エポック数
	constexpr size_t EpochDisp = 10;	// エポック表示数
	constexpr ActvFunc InFunc  = ActvFunc::ReLU;					// 入力層の活性化関数
	constexpr ActvFunc HidFunc = ActvFunc::ReLU;					// 隠れ層の活性化関数
	constexpr ActvFunc OutFunc = ActvFunc::SOFTMAX;					// 出力層の活性化関数
	constexpr NnInitTypes InitType = NnInitTypes::HE;				// 重み初期化のタイプ
	constexpr NnDescentTypes GradDesType = NnDescentTypes::MOMENTUM;// 勾配降下法のタイプ
	constexpr NnDropout DropoutEna = NnDropout::DISABLE;			// ドロップアウトするかどうか
	constexpr NnShuffle EpochShfl = NnShuffle::DISABLE;				// エポック毎にデータセットをシャッフルするかどうか
	constexpr double epsilon = 0.001;	// 学習率
	constexpr double alpha = 0.9;		// 学習率
	//constexpr double beta = 0.999;		// 学習率(未使用)
	//constexpr double zero = 1e-8;		// ゼロ割回避用の係数(未使用)
	constexpr double p = 0.95;			// ドロップアウト率
	
	// 3層ニューラルネットワークの生成
	FeedforwardNeuralNet3<
		IrisClassDatasets<MinbatSize>,	// データセットクラス名
		InSize,			// 入力層のサイズ
		HidSize,		// 隠れ層のサイズ
		OutSize,		// 出力層のサイズ 
		MinbatSize,		// ミニバッチサイズ
		MinbatNum,		// ミニバッチ数
		Epoch,			// エポック数
		EpochDisp,		// エポック表示数
		InFunc,			// 入力層の活性化関数
		HidFunc,		// 隠れ層の活性化関数
		OutFunc,		// 出力層の活性化関数
		InitType,		// 重み初期化のタイプ
		GradDesType,	// 勾配降下法のタイプ
		DropoutEna,		// ドロップアウトイネーブル
		EpochShfl		// エポック毎データセットシャッフル
	> FFNN;
	
	// フィードフォワードニューラルネットワークの訓練
	IrisClassDatasets<MinbatSize> IrisData;		// アヤメデータセットの生成(分類確率ベクトル版)
	IrisData.DispMeasuredData();				// 計測データの表示
	IrisData.DispClassData();					// 分類データの表示
	FFNN.SetGainOfMomentumSGD(epsilon, alpha);	// モーメンタム勾配降下法の学習率の設定
	FFNN.SetDropoutRate(p);						// ドロップアウトの設定
	FFNN.Train(IrisData);						// ニューラルネットの訓練
	FFNN.DispWeight();							// 学習後の重み行列の表示
	FFNN.DispFinalTestData();					// 順伝播出力値と正解値の表示
	FFNN.WriteLeaningCurvePNG("LearningCurve.png", -50, 10);	// 学習曲線のPNGファイル出力
	
	// パーセプトロンの設定値，重み行列，バイアスベクトルをCSVファイルに保存
	FFNN.SaveSettings("FFNN");
	FFNN.SaveWeightAndBias("FFNN");
	
	// CSVファイルからパーセプトロンの設定値，重み行列，バイアスベクトルに読み込む関数
	FFNN.LoadSettings("FFNN");
	FFNN.LoadWeightAndBias("FFNN");
	
	// パーセプトロン設定値，重み行列，バイアスベクトルの表示
	FFNN.DispSettings();
	FFNN.DispWeight();
	FFNN.DispBias();
	
	// 訓練済みニューラルネットワークを使った推定
	Matrix<1,InSize> x;		// ニューラルネットワークの入力ベクトル
	Matrix<1,OutSize> y;	// ニューラルネットワークの出力ベクトル
	x.Set(5.9,3.0,5.1,1.8);	// ニューラルネットワークへの入力値(アヤメの計測データ)
	FFNN.Estimate(x, y);	// ニューラルネットワークの推定計算(アヤメの種類の推測)
	PrintMatrix(y, "% g");	// ニューラルネットワークの出力の表示([0,0,1]^Tが表示されればバージニカなので正解)
	
	return EXIT_SUCCESS;	// 正常終了
}

