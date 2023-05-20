//! @file FeedforwardNeuralNet3.hh
//! @brief 順伝播型ニューラルネットワーク（3層版）
//!
//! 順伝播型ニューラルネットワーク（3層版）
//!
//! @date 2020/05/17
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef FEEDFORWARDNEURALNET3
#define FEEDFORWARDNEURALNET3

#include <cassert>
#include <array>
#include <string>
#include "Matrix.hh"
#include "NeuralNetParamDef.hh"
#include "SingleLayerPerceptron.hh"
#include "FrameGraphics.hh"
#include "CuiPlot.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

namespace ARCS {	// ARCS名前空間
	//! @brief 多層パーセプトロンクラス3層版
	// 下記のテンプレートの実装はとてもアホっぽいが，
	// C++17ではテンプレートパラメータにローカルなstructも配列も渡せないし，
	// パラメータパックは展開めんどくせーし，という言い訳をここに書いておく
	//! @tparam	DsName	データセットクラス名
	//! @tparam	EpochDisp	エポックループ表示数
	//! @tparam	DsShfl	エポック毎にシャッフルするかどうかの設定
	template <
		typename DsName,	// データセットクラス名
		size_t InSize,		// 入力層のサイズ
		size_t HidSize,		// 隠れ層のサイズ
		size_t OutSize,		// 出力層のサイズ 
		size_t MinbatSize,	// ミニバッチサイズ
		size_t MinbatNum,	// ミニバッチ数
		size_t Epoch,		// エポック数
		size_t EpochDisp,	// エポックループ表示数
		ActvFunc InFunc,			// 入力層の活性化関数
		ActvFunc HidFunc,			// 隠れ層の活性化関数
		ActvFunc OutFunc,			// 出力層の活性化関数
		NnInitTypes InitType,		// 重み初期化のタイプ
		NnDescentTypes GradDesType,	// 勾配降下法のタイプ
		NnDropout DropOutEnable,	// ドロップアウトイネーブル
		NnShuffle DsShfl			// エポック毎のデータセットシャッフル
	>
	class FeedforwardNeuralNet3 {
		public:
			//! @brief コンストラクタ
			FeedforwardNeuralNet3()
				: InputLayer(), HiddenLayer(), OutputLayer(),
				  X(), ZIn(), ZHid(), Y(), D(), WDeltaOut(), WDeltaHid(), WDeltaIn(),
				  EpochNumbers(), TrainLoss(), TestLoss(),
				  XTest(), ZInTest(), ZHidTest(), YTest(), DTest(),
				  zin(), zhid()
			{
				PassedLog();
			}
			
			//! @brief ムーブコンストラクタ
			//! @param[in]	r	右辺値
			FeedforwardNeuralNet3(FeedforwardNeuralNet3&& r)
				: InputLayer(r.InputLayer), HiddenLayer(r.HiddenLayer), OutputLayer(r.OutputLayer),
				  X(r.X), ZIn(r.ZIn), ZHid(r.ZHid), Y(r.Y), D(r.D), WDeltaOut(r.WDeltaOut), WDeltaHid(r.WDeltaHid), WDeltaIn(r.WDeltaIn),
				  EpochNumbers(r.EpochNumbers), TrainLoss(r.TrainLoss), TestLoss(r.TestLoss),
				  XTest(r.XTest), ZInTest(r.ZInTest), ZHidTest(r.ZHidTest), YTest(r.YTest), DTest(r.DTest),
				  zin(r.zin), zhid(r.zhid)
			{
				
			}
			
			//! @brief デストラクタ
			~FeedforwardNeuralNet3(){
				PassedLog();
			}
			
			//! @brief 重み行列の初期化
			void InitWeight(void){
				InputLayer.InitWeight(InSize);	// 入力層の重み行列の乱数による初期化
				HiddenLayer.InitWeight(InSize);	// 内部層の重み行列の乱数による初期化
				OutputLayer.InitWeight(HidSize);// 出力層の重み行列の乱数による初期化
			}
			
			//! @brief モーメンタム確率的勾配降下法の更新ゲイン(学習率)の設定
			//! @param[in]	epsilon	更新ゲイン(学習率)
			//! @param[in]	alpha	更新ゲイン(運動量項の学習率)
			void SetGainOfMomentumSGD(double epsilon, double alpha){
				InputLayer.SetGainOfMomentumSGD(epsilon, alpha);
				HiddenLayer.SetGainOfMomentumSGD(epsilon, alpha);
				OutputLayer.SetGainOfMomentumSGD(epsilon, alpha);
			}
			
			//! @brief ドロップアウト率の設定
			//! @param[in]	DropoutRate	ドロップアウト率(1のときドロップアウトしない，0のとき全部ドロップアウト)
			void SetDropoutRate(double DropoutRate){
				InputLayer.SetDropoutRate(DropoutRate);
				HiddenLayer.SetDropoutRate(DropoutRate);
				OutputLayer.SetDropoutRate(DropoutRate);
			}
			
			//! @brief 誤差逆伝播法による訓練をする関数
			//! @param[in]	Datasets	データセット
			void Train(DsName& Datasets){
				InitWeight();				// 重み行列の初期化
				InputLayer.NomalizeDataset(Datasets.MeasuredData);				// データセットの正規化(標準化)
				Datasets.GetMeasuredBatchData(Datasets.FinalMinbatNum, XTest);	// 最後のミニバッチをテスト入力データとして取得
				Datasets.GetClassBatchData(Datasets.FinalMinbatNum, DTest);		// 最後のミニバッチをテスト目標データとして取得
				
				// 誤差逆伝播法のエポック数分のループ
				printf("\nBackpropagation Training:\n");
				printf("Epoch : Train Loss [dB]    Test Loss [dB]\n");
				for(size_t i = 0; i < Epoch; ++i){
					// ミニバッチ数分のループ
					for(size_t j = 0; j < MinbatNum; ++j){
						// 訓練データの準備
						Datasets.GetMeasuredBatchData(j + 1, X);// ミニバッチ分の入力データを取得
						Datasets.GetClassBatchData(j + 1, D);	// ミニバッチ分の目標データを取得
						
						// ドロップアウトの準備
						InputLayer.CalcDropout();
						HiddenLayer.CalcDropout();
						
						// 順伝播計算
						InputLayer.CalcForwardForTraining(X, ZIn);		// 入力層の順伝播計算
						HiddenLayer.CalcForwardForTraining(ZIn, ZHid);	// 内部層の順伝播計算
						OutputLayer.CalcForwardForTraining(ZHid, Y);	// 出力層の順伝播計算
						
						// 逆伝播による誤差行列計算
						OutputLayer.CalcDeltaForOutputLayer(Y, D, WDeltaOut);	// 出力層の逆伝播計算
						HiddenLayer.CalcDelta(WDeltaOut, WDeltaHid);			// 内部層の逆伝播計算
						InputLayer.CalcDelta(WDeltaHid, WDeltaIn);				// 入力層の逆伝播計算
						
						// 重み行列の更新
						InputLayer.UpdateWeight(X);		// 入力層の重み更新
						HiddenLayer.UpdateWeight(ZIn);	// 内部層の重み更新
						OutputLayer.UpdateWeight(ZHid);	// 出力層の重み更新
					}
					if constexpr(DsShfl == NnShuffle::ENABLE) Datasets.ShuffleDatasets();	// データセットのシャッフル
					
					// テストデータの順伝播計算
					InputLayer.CalcForwardForEstimation(XTest, ZInTest);	// 入力層の順伝播計算
					HiddenLayer.CalcForwardForEstimation(ZInTest, ZHidTest);// 内部層の順伝播計算
					OutputLayer.CalcForwardForEstimation(ZHidTest, YTest);	// 出力層の順伝播計算
					
					// 学習曲線評価用のデータ保持
					EpochNumbers[i] = i;
					TrainLoss[i] = 10.0*log10(OutputLayer.GetLoss(Y, D));		// 訓練誤差
					TestLoss[i] = 10.0*log10(OutputLayer.GetLoss(YTest, DTest));// テスト誤差
					
					// エポック毎の訓練誤差とテスト誤差の表示
					if(i % (Epoch/EpochDisp) == 0){
						printf("%5lu : % 16.8f, % 16.8f\n", i, TrainLoss[i], TestLoss[i]);
					}
				}
			}
			
			//! @brief 重み行列の表示
			void DispWeight(void){
				printf("\nWeight Matrices of All Layers:\n");
				InputLayer.DispWeight();	// 入力層の重みの表示
				HiddenLayer.DispWeight();	// 内部層の重みの表示
				OutputLayer.DispWeight();	// 出力層の重みの表示
			}
			
			//! @brief バイアスベクトルの表示
			void DispBias(void){
				printf("\nBias Vectors of All Layers:\n");
				InputLayer.DispBias();
				HiddenLayer.DispBias();
				OutputLayer.DispBias();
			}
			
			//! @brief パーセプトロン設定値の表示
			void DispSettings(void){
				printf("\nNeural Network Settings:\n");
				InputLayer.DispSettings();
				HiddenLayer.DispSettings();
				OutputLayer.DispSettings();
			}
			
			//! @brief 最後の学習結果を使った順伝播出力値と正解値の確認
			void DispFinalTestData(void){
				printf("\nFinal Test Data Comfirmation:\n");
				PrintMatrix(YTest, "% g");	// 順伝播出力値
				PrintMatrix(DTest, "% g");	// 正解値
			}
			
			//! @brief 学習曲線グラフのPNG画像ファイルを書き出す関数
			//! @param[in]	FileName	ファイル名
			//! @param[in]	YMIN	グラフの縦軸最小値 [dB]
			//! @param[in]	YMAX	グラフの縦軸最大値 [dB]
			void WriteLeaningCurvePNG(const std::string& FileName, const int Ymin, const int Ymax){
				// 学習曲線のグラフの表示
				FrameGraphics FG(GRAPH_WIDTH, GRAPH_HEIGHT);
				CuiPlot Plot(FG, 0, 0, GRAPH_WIDTH, GRAPH_HEIGHT);
				Plot.SetAxisLabels("Epoch Number", "Errors 10log10 [dB]");
				Plot.SetRanges(0, Epoch, (double)Ymin, (double)Ymax);
				Plot.SetGridLabelFormat("%5.0f", "%2.1f");
				Plot.DrawAxis();
				Plot.DrawLegend(1, "Training", FGcolors::CYAN);
				Plot.DrawLegend(2, "Test", FGcolors::MAGENTA);
				Plot.Plot(EpochNumbers, TrainLoss, CuiPlotTypes::PLOT_BOLDSTAIRS, FGcolors::CYAN);
				Plot.Plot(EpochNumbers, TestLoss, CuiPlotTypes::PLOT_BOLDSTAIRS, FGcolors::MAGENTA);
				FG.SavePngImageFile(FileName);
			}
			
			//! @brief 訓練済みニューラルネットワークを使った推定計算
			//! @param[in]	x	ニューラルネットワークの入力ベクトル
			//! @param[out]	y	ニューラルネットワークの出力ベクトル
			void Estimate(const Matrix<1,InSize>& x, Matrix<1,OutSize>& y){
				// 順伝播計算
				Matrix<1,InSize> xn = x;			// 正規化用
				InputLayer.NormalizeInput(xn);		// 入力データの正規化(標準化)
				InputLayer.CalcForwardForEstimation(xn, zin);	// 入力層の順伝播計算
				HiddenLayer.CalcForwardForEstimation(zin, zhid);// 内部層の順伝播計算
				OutputLayer.CalcForwardForEstimation(zhid, y);	// 出力層の順伝播計算
			}
			
			//! @brief 各レイヤの重み行列とバイアスベクトルをCSVファイルとして保存する関数
			//! @param[in]	FileName	CSVファイル名(拡張子なし)
			void SaveWeightAndBias(const std::string& FileName){
				InputLayer.SaveWeightAndBias(FileName + "-InputLayerW.csv", FileName + "-InputLayerb.csv");
				HiddenLayer.SaveWeightAndBias(FileName + "-HiddenLayerW.csv", FileName + "-HiddenLayerb.csv");
				OutputLayer.SaveWeightAndBias(FileName + "-OutputLayerW.csv", FileName + "-OutputLayerb.csv");
			}
			
			//! @brief CSVファイルから各レイヤの重み行列とバイアスベクトルに読み込む関数
			//! @param[in]	FileName	CSVファイル名(拡張子なし)
			void LoadWeightAndBias(const std::string& FileName){
				InputLayer.LoadWeightAndBias(FileName + "-InputLayerW.csv", FileName + "-InputLayerb.csv");
				HiddenLayer.LoadWeightAndBias(FileName + "-HiddenLayerW.csv", FileName + "-HiddenLayerb.csv");
				OutputLayer.LoadWeightAndBias(FileName + "-OutputLayerW.csv", FileName + "-OutputLayerb.csv");
			}
			
			//! @brief 各レイヤのパーセプトロンの設定値をCSVファイルとして保存する関数
			//! @param[in]	FileName	CSVファイル名(拡張子なし)
			void SaveSettings(const std::string& FileName){
				InputLayer.SaveSettings(FileName + "-InputLayerSet.csv");
				HiddenLayer.SaveSettings(FileName + "-HiddenLayerSet.csv");
				OutputLayer.SaveSettings(FileName + "-OutputLayerSet.csv");
			}
			
			//! @brief CSVファイルから各レイヤのパーセプトロンの設定値を読み込む関数
			//! @param[in]	FileName	CSVファイル名(拡張子なし)
			void LoadSettings(const std::string& FileName){
				InputLayer.LoadSettings(FileName + "-InputLayerSet.csv");
				HiddenLayer.LoadSettings(FileName + "-HiddenLayerSet.csv");
				OutputLayer.LoadSettings(FileName + "-OutputLayerSet.csv");
			}
			
		private:
			FeedforwardNeuralNet3(const FeedforwardNeuralNet3&) = delete;					//!< コピーコンストラクタ使用禁止
			const FeedforwardNeuralNet3& operator=(const FeedforwardNeuralNet3&) = delete;	//!< 代入演算子使用禁止
			
			// 学習曲線グラフの設定
			static constexpr int GRAPH_WIDTH = 1000;	//!< [px] グラフの横幅
			static constexpr int GRAPH_HEIGHT = 500;	//!< [px] グラフの高さ
			
			// 3層ニューラルネットワーク
			SingleLayerPerceptron<InSize,  InSize,  MinbatSize, InFunc,  InitType, GradDesType, DropOutEnable> InputLayer;	//!< 入力層
			SingleLayerPerceptron<InSize,  HidSize, MinbatSize, HidFunc, InitType, GradDesType, DropOutEnable> HiddenLayer;	//!< 隠れ層
			SingleLayerPerceptron<HidSize, OutSize, MinbatSize, OutFunc, InitType, GradDesType> OutputLayer;				//!< 出力層
			
			// 誤差逆伝播用の変数
			Matrix<MinbatSize, InSize> X;			//!< 順伝播入力行列(ミニバッチ)
			Matrix<MinbatSize, InSize> ZIn;			//!< 入力層の出力行列
			Matrix<MinbatSize, HidSize> ZHid;		//!< 内部層の出力行列
			Matrix<MinbatSize, OutSize> Y;			//!< 順伝播出力行列
			Matrix<MinbatSize, OutSize> D;			//!< 目標値行列
			Matrix<MinbatSize, HidSize> WDeltaOut;	//!< 出力層の重み誤差行列
			Matrix<MinbatSize, InSize> WDeltaHid;	//!< 内部層の重み誤差行列
			Matrix<MinbatSize, InSize> WDeltaIn;	//!< 入力層の重み誤差行列
			std::array<double, Epoch> EpochNumbers;	//!< エポック番号
			std::array<double, Epoch> TrainLoss;	//!< 訓練誤差
			std::array<double, Epoch> TestLoss;		//!< テスト誤差
			
			// テストデータ用の変数
			Matrix<MinbatSize, InSize> XTest;		//!< テスト用順伝播入力行列(ミニバッチ)
			Matrix<MinbatSize, InSize> ZInTest;		//!< 入力層の出力行列
			Matrix<MinbatSize, HidSize> ZHidTest;	//!< 内部層の出力行列
			Matrix<MinbatSize, OutSize> YTest;		//!< 順伝播出力行列
			Matrix<MinbatSize, OutSize> DTest;		//!< 目標値行列
			
			// 推定計算用の変数
			Matrix<1,InSize> zin;	//!< 入力層の出力ベクトル
			Matrix<1,HidSize> zhid;	//!< 内部層の出力ベクトル
	};

}

#endif

