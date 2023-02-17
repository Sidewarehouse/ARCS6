//! @file CurrencyDatasets.cc
//! @brief 機械学習用 為替データセットクラス
//!
//! 機械学習のための為替データを扱うためのデータセットクラス
//! データ自体は保持していないので，外部のCSVデータを必要とする
//! 注意：CSVデータの改行コードは「LF」とし、「CRLF」は不可
//!
//! @date 2020/05/25
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef CURRENCYDATASETS
#define CURRENCYDATASETS

#include <cassert>
#include <array>
#include <string>
#include <ctime>
#include "Matrix.hh"
#include "BatchProcessor.hh"
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
//! @brief 機械学習用 為替データセットクラス
//! @tparam	D	データ数
//! @tparam	M	ミニバッチ数
template <size_t D, size_t M>
class CurrencyDatasets {
	public:
		// データセット情報の定義
		static constexpr size_t N = 1;	//!< 計測データチャネル数
		static constexpr size_t C = 1;	//!< 訓練データチャネル数
		static constexpr size_t FinalMinbatNum = D/M - 1;	//!< 最後のミニバッチ番号
		
		// 為替データ
		Matrix<N,D> CurrencyData;				//!< 為替データ
		Matrix<N,D> TimeData;					//!< Unix時間ベクトル
		std::array<std::string, D> TimeString;	//!< 時刻文字データ
		
		//! @brief コンストラクタ
		//! @param[in]	CurrencyFileName	為替データのCSVファイル名
		//! @param[in]	TimedataFileName	時間データのCSVファイル名
		CurrencyDatasets(const std::string& CurrencyFileName, const std::string& TimedataFileName)
			: CurrencyData(), TimeData(), TimeString()
		{
			PassedLog();
			CsvManipulator::LoadFile(CurrencyData, CurrencyFileName);	// 為替データを行列として読み込み
			CsvManipulator::LoadFile(TimeString, TimedataFileName);		// 時間データを文字列配列として読み込み
			//ConvTimeStringToData();		// 時刻文字データをUnix時間に変換
			TimeData = Matrix<1,D>::ramp();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		CurrencyDatasets(CurrencyDatasets&& r)
			: CurrencyData(r.CurrencyData), TimeData(r.TimeData), TimeString(r.TimeString)
		{
			
		}
		
		//! @brief デストラクタ
		~CurrencyDatasets(){
			PassedLog();
		}
		
		//! @brief 為替データを表示する関数
		//! @param[in]	DispNum	表示間引き数
		void DispCurrencyData(const size_t DispNum){
			printf("\nCurrency Data:\n");
			for(size_t i = 0; i < D; ++i){
				if((i % DispNum) == 0) printf("%s : %10.0f : %7.3f\n", TimeString.at(i).c_str(), TimeData[i+1], CurrencyData[i+1]);
			}
		}
		
		//! @brief 為替データの時系列プロットを表示する関数
		void WritePngCurrencyPlot(const double Min, const double Max, const std::string& FileName){
			FrameGraphics FG(GRAPH_WIDTH, GRAPH_HEIGHT);
			CuiPlot Plot(FG, 0, 0, GRAPH_WIDTH, GRAPH_HEIGHT);
			Plot.SetAxisLabels("Time Index", FileName);
			Plot.SetRanges(0, D, Min, Max);
			Plot.SetGridLabelFormat("%5.0f", "%3.0f");
			Plot.DrawAxis();
			Plot.DrawLegend(1, "USDJPY", FGcolors::CYAN);
			Plot.Plot(TimeData, CurrencyData, CuiPlotTypes::PLOT_LINE, FGcolors::CYAN);
			FG.SavePngImageFile(FileName);
		}
		
		/*
		//! @brief i番目のミニバッチ為替データを取得する関数
		//! @param[in]	i	ミニバッチ番号
		//! @param[out]	Y	為替データ
		void GetCurrectBatchData(const size_t i, Matrix<M,N>& Y){
			arcs_assert(i <= FinalMinbatNum);	// ミニバッチ番号の範囲チェック
			Y = BatchProcessor::GetMiniBatchData<N,D,M>(MeasuredData, i);	// データセットを切り出してミニバッチにする
		}
		*/
		
	private:
		CurrencyDatasets(const CurrencyDatasets&) = delete;					//!< コピーコンストラクタ使用禁止
		const CurrencyDatasets& operator=(const CurrencyDatasets&) = delete;//!< 代入演算子使用禁止
		
		// 学習曲線グラフの設定
		static constexpr int GRAPH_WIDTH = 1000;	//!< [px] グラフの横幅
		static constexpr int GRAPH_HEIGHT = 500;	//!< [px] グラフの高さ
		
		//! @brief 時刻文字データからUnix時間データへ変換する関数 (遅い！)
		void ConvTimeStringToData(void){
			struct tm TimeBuff;		// 時刻構造体バッファ
			time_t UnixTime;		// Unix時間
			
			// 1回目がなぜかおかしくなるので空回ししておく
			strptime(TimeString.at(0).c_str(), "%Y.%m.%d.%H:%M", &TimeBuff);	// 文字列をパースして時刻構造体に変換
			UnixTime = mktime(&TimeBuff);
			
			// 全時刻データをUnix時間に変換
			for(size_t i = 0; i < D; ++i){
				strptime(TimeString.at(i).c_str(), "%Y.%m.%d.%H:%M", &TimeBuff);	// 文字列をパースして時刻構造体に変換
				UnixTime = mktime(&TimeBuff);
				TimeData[i+1] = (double)UnixTime;
			}
		}
};
}

#endif

