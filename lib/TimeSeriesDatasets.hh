//! @file TimeSeriesDatasets.hh
//! @brief 機械学習用 時系列データセットクラス
//!
//! 機械学習のための時系列データを扱うためのデータセットクラス
//! データ自体は保持していないので，外部のCSVデータを必要とする
//! 注意：CSVデータの改行コードは「LF」とし、「CRLF」は不可
//!
//! @date 2021/08/19
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TIMESERIESDATASETS
#define TIMESERIESDATASETS

#include <cassert>
#include <array>
#include <string>
#include <ctime>
#include "Matrix.hh"
#include "Statistics.hh"
#include "CsvManipulator.hh"
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
//! @brief 機械学習用 時系列データセットクラス
//! @tparam N	入力データチャネル数
//! @tparam K	訓練データチャネル数
//! @tparam	T	時間方向のデータ数
//! @tparam	W	入力ウィンドウ幅
//! @tparam	M	ミニバッチサイズ
template <size_t N, size_t K, size_t T, size_t W, size_t M>
class TimeSeriesDatasets {
	public:
		// 時系列データ
		std::array<Matrix<M,N>, W + 2> InputData;//!< ベクトル配列版の標準化済み入力データ(範囲 t = 1 … W, t = 0 と W + 1 の分も確保)
		Matrix<M,K> TrainData;//!< ベクトル配列版の標準化済み訓練データ(範囲 t = 1 … W, t = 0 と W + 1 の分も確保)
		
		//! @brief 空コンストラクタ
		TimeSeriesDatasets(void)
			: InputData(), TrainData(), 
			  TimeStamp(), InputMat(), TrainMat(), StdInputMat(), StdTrainMat(),
			  xbarIn(), sigmaIn(), xbarTr(), sigmaTr()
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ
		//! @param[in]	InputFile	時系列入力データのCSVファイル名
		//! @param[in]	TrainFile	時系列訓練データのCSVファイル名
		TimeSeriesDatasets(const std::string& InputFile, const std::string& TrainFile)
			: InputData(), TrainData(), 
			  TimeStamp(), InputMat(), TrainMat(), StdInputMat(), StdTrainMat(),
			  xbarIn(), sigmaIn(), xbarTr(), sigmaTr()
		{
			PassedLog();
			TimeStamp = tp(Matrix<1,T>::ramp());						// タイムスタンプを生成
			
			// 時系列入力データの前処理
			Matrix<N,T> InputDataBuf;									// 入力データバッファ
			CsvManipulator::LoadFile(InputDataBuf, InputFile);			// CSVデータを行列として読み込み
			InputMat = tp(InputDataBuf);								// 今後のために転置をしておく
			StandardizeDataset(InputMat, xbarIn, sigmaIn, StdInputMat);	// データセットの標準化
			VectorizeInputData(StdInputMat, InputData);					// 縦ベクトル配列化
			
			// 時系列訓練データの前処理
			Matrix<K,T> TrainDataBuf;									// 訓練データバッファ
			CsvManipulator::LoadFile(TrainDataBuf, TrainFile);			// CSVデータを行列として読み込み
			TrainMat = tp(TrainDataBuf);								// 今後のために転置をしておく
			StandardizeDataset(TrainMat, xbarTr, sigmaTr, StdTrainMat);	// データセットの標準化
			VectorizeTrainData(StdTrainMat, TrainData);					// 縦ベクトル配列化
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		TimeSeriesDatasets(TimeSeriesDatasets&& r)
			: TimeStamp(r.TimeStamp), InputMat(r.InputTimeData), TrainMat(r.TrainMat),
			  StdInputMat(r.StdInputMat), StdTrainMat(r.StdTrainMat),
			  InputData(r.InputData), TrainData(r.TrainData),
			  xbarIn(r.xbarIn), sigmaIn(r.sigmaIn), xbarTr(r.xbarTr), sigmaTr(r.sigmaTr)
		{
			
		}
		
		//! @brief デストラクタ
		~TimeSeriesDatasets(){
			PassedLog();
		}
		
		//! @brief 入力データの標準化
		//! @param[in]	x	入力データ
		void StandardizeInput(Matrix<1,N>& x){
			x = (x - xbarIn) % sigmaIn;	// 入力データの標準化値に基づいて入力を標準化
		}
		
		//! @brief 標準化済み入力＆訓練データを表示する関数
		//! @param[in]	DispNum	表示間引き数
		void DispInputAndTrainData(const size_t DispNum){
			printf("\nTime Series Input And Train Data :\n");
			for(size_t i = 1; i <= W; ++i){
				//if((i % DispNum) == 0){
					printf("t = %zu :\n", i);
					PrintMat(InputData.at(i));
				//}
			}
			PrintMat(TrainData);
			printf("Standardization Info :\n");
			PrintMat(xbarIn);
			PrintMat(sigmaIn);
			PrintMat(xbarTr);
			PrintMat(sigmaTr);
		}
		
		//! @brief 標準化済み入力データの時間波形をPNGファイルに書き出す関数
		//! @param[in]	Min	グラフの最小値
		//! @param[in]	Max	グラフの最大値
		//! @param[in]	FileName	出力するファイル名
		void WriteInputPlot(const double Min, const double Max, const std::string& FileName){
			// グラフ設定
			FrameGraphics FG(GRAPH_WIDTH, GRAPH_HEIGHT);
			CuiPlot Plot(FG, 0, 0, GRAPH_WIDTH, GRAPH_HEIGHT);
			Plot.SetAxisLabels("Time Index", FileName);
			Plot.SetRanges(0, T, Min, Max);
			Plot.SetGridLabelFormat("%5.0f", "%3.0f");
			Plot.DrawAxis();
			
			// 入力データチャネルごとにプロット
			for(size_t i = 1; i <= N; ++ i){
				Plot.DrawLegend(i, "Variable_i", FGcolors::CYAN);
				Plot.Plot(tp(TimeStamp), tp(getrow(StdInputMat,i)), CuiPlotTypes::PLOT_LINE, FGcolors::CYAN);
			}
			
			// PNG画像ファイル出力
			FG.SavePngImageFile(FileName);
		}
		
		//! @brief 標準化済み入力データの時間波形をPNGファイルに書き出す関数
		//! @param[in]	Min	グラフの最小値
		//! @param[in]	Max	グラフの最大値
		//! @param[in]	FileName	出力するファイル名
		void WriteTrainPlot(const double Min, const double Max, const std::string& FileName){
			// グラフ設定
			FrameGraphics FG(GRAPH_WIDTH, GRAPH_HEIGHT);
			CuiPlot Plot(FG, 0, 0, GRAPH_WIDTH, GRAPH_HEIGHT);
			Plot.SetAxisLabels("Time Index", FileName);
			Plot.SetRanges(0, T, Min, Max);
			Plot.SetGridLabelFormat("%5.0f", "%3.0f");
			Plot.DrawAxis();
			
			// 入力データチャネルごとにプロット
			for(size_t i = 1; i <= K; ++ i){
				Plot.DrawLegend(i, "Variable_i", FGcolors::CYAN);
				Plot.Plot(tp(TimeStamp), tp(getrow(StdTrainMat,i)), CuiPlotTypes::PLOT_LINE, FGcolors::MAGENTA);
			}
			
			// PNG画像ファイル出力
			FG.SavePngImageFile(FileName);
		}
		
	private:
		TimeSeriesDatasets(const TimeSeriesDatasets&) = delete;					//!< コピーコンストラクタ使用禁止
		const TimeSeriesDatasets& operator=(const TimeSeriesDatasets&) = delete;//!< 代入演算子使用禁止
		
		// 時系列データの時間波形グラフの設定
		static constexpr int GRAPH_WIDTH = 1000;	//!< [px] グラフの横幅
		static constexpr int GRAPH_HEIGHT = 500;	//!< [px] グラフの高さ
		
		// 内部用時系列データ
		Matrix<T,1> TimeStamp;	//!< タイムスタンプ
		Matrix<T,N> InputMat;	//!< 時系列入力データ(横方向を時刻，縦方向をチャネル数とする)
		Matrix<T,K> TrainMat;	//!< 時系列訓練データ(横方向を時刻，縦方向をチャネル数とする)
		Matrix<T,N> StdInputMat;//!< 行列版の標準化済み入力データ
		Matrix<T,K> StdTrainMat;//!< 行列版の標準化済み訓練データ
		
		// データセットの各種パラメータ
		Matrix<1,N> xbarIn;		//!< 入力データ標準化用の平均値
		Matrix<1,N> sigmaIn;	//!< 入力データ標準化用の標準偏差
		Matrix<1,N> xbarTr;		//!< 訓練データ標準化用の平均値
		Matrix<1,N> sigmaTr;	//!< 訓練データ標準化用の標準偏差
		
		//! @brief 生計測データセットを標準化する関数
		//! @param[in]	RawData	生データ
		//! @param[in]	Mean	横方向の平均値(縦ベクトル)
		//! @param[in]	Std		横方向の標準偏差(縦ベクトル)
		static void StandardizeDataset(const Matrix<T,N>& RawData, Matrix<1,N>& Mean, Matrix<1,N>& Std, Matrix<T,N>& StdData){
			Statistics::MeanRow(RawData, Mean);				// 行方向の平均の計算
			Statistics::StandardDeviationRow(RawData, Std);	// 行方向の標準偏差の計算
			Matrix<1,N> buff;
			for(size_t i = 1; i <= T; ++i){
				buff = (getcolumn(RawData, i) - Mean) % Std;// データセットを標準化
				setcolumn(StdData, buff, i);
			}
		}
		
		//! @brief	入力データを縦ベクトル配列化する関数
		//! @param[in]	MatData	入力行列データ
		//! @param[out]	VecData	配列化された縦ベクトル入力データ
		static void VectorizeInputData(const Matrix<T,N>& MatData, std::array<Matrix<M,N>, W + 2>& VecData){
			for(size_t i = 1; i <= W; ++i){
				getvvector(MatData, i, 1, VecData.at(i));
			}
		}
		
		//! @brief	訓練データを縦ベクトル配列化する関数
		//! @param[in]	MatData	訓練行列データ
		//! @param[out]	VecData	配列化された縦ベクトル訓練データ
		static void VectorizeTrainData(const Matrix<T,K>& MatData, Matrix<M,K>& VecData){
			getvvector(MatData, W + 1, 1, VecData);
		}
		
};
}

#endif

