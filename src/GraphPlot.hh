//! @file GraphPlot.cc
//! @brief グラフプロットクラス
//!
//! グラフを描画するクラス
//!
//! @date 2023/05/20
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef GRAPHPLOT
#define GRAPHPLOT

#include <pthread.h>
#include <cfloat>
#include "ConstParams.hh"
#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
//! @brief グラフプロットクラス
class GraphPlot {
	public:
		CuiPlot UserPlot;			//!< ユーザプロット
		
		explicit GraphPlot(void);	//!< コンストラクタ
		~GraphPlot();				//!< デストラクタ
		void DrawPlotPlane(void);	//!< プロット平面の描画
		void DrawWaves(void);		//!< プロット波形の描画
		void ResetWaves(void);		//!< プロットをリセットする関数
		void SaveScreenImage(void);	//!< 画面をPNGファイルとして出力する関数
		
		//! @brief プロット描画時間に値を設定する関数
		//! @param[in]	T	周期 [s]
		//! @param[in]	t	時刻 [s]
		void SetTime(const double T, const double t){
			const double tlimited = fmod(t, ConstParams::PLOT_TIMESPAN);	// [s] 横軸時間の範囲を[0～最大の時刻]に留める計算
			const double tstorage = fmod(t, ConstParams::PLOT_TIMERESO);	// [s] リングバッファ保存時間になったかの判定用時刻の計算
			if(tstorage <= T){
				// リングバッファ保存時刻になったら
				StorageEnable = true;					// リングバッファ保存を有効にして
				TimeRingBuf.SetFirstValue(tlimited);	// 時刻をリングバッファに詰める
			}else{
				// リングバッファ保存時刻以外の場合は
				StorageEnable = false;	// 保存を無効にする
			}
		}
		
		//! @brief プロット描画変数に値を設定する関数(可変長引数テンプレート)
		//! @param[in] u1...u2 インジケータの値
		template<typename T1, typename... T2>
		void SetVars(const T1& u1, const T2&... u2){
			// リングバッファ保存時刻ではなかったら
			if(StorageEnable == false) return;	// 何もせずに終了
			
			// 再帰で順番に可変長引数を読み込んでいく
			if(VarsCount == 0){
				PlotNumBuf = (size_t)u1;	// 1個目の引数はプロット番号として格納
			}else{
				if(VarsCount <= ConstParams::PLOT_VAR_NUM[PlotNumBuf]){
					// 変数要素数が有効な範囲内なら
					//pthread_mutex_lock(&PlotVarsMutex); ←リアルタイム性悪化の原因！(一時的な対処)
					VarsRingBuf.at(PlotNumBuf).at(VarsCount - 1).SetFirstValue((double)u1);	// 変数値リングバッファに詰める
					//pthread_mutex_unlock(&PlotVarsMutex);
				}
			}
			
			++VarsCount;	// 再帰カウンタをインクリメント
			SetVars(u2...);	// 自分自身を呼び出す(再帰)
		}
		//! @brief 再帰の最後に呼ばれる関数
		void SetVars(){
			VarsCount = 0;	// すべての作業が終わったので，再帰カウンタを零に戻しておく
		}
		
		//! @brief 作業空間プロットに位置ベクトルを設定する関数（6軸配列版）
		void SetWorkspace(const std::array<Matrix<1,6>, 6>& AxPosition);
		//! @brief 作業空間プロットに位置ベクトルを設定する関数（3軸個別版）
		void SetWorkspace(const Matrix<1,3>& Pos1, const Matrix<1,3>& Pos2, const Matrix<1,3>& Pos3);
		
	private:
		GraphPlot(const GraphPlot&) = delete;					//!< コピーコンストラクタ使用禁止
		const GraphPlot& operator=(const GraphPlot&) = delete;	//!< 代入演算子使用禁止
		GraphPlot(GraphPlot&& r) = delete;						//!< ムーブコンストラクタ使用禁止
		
		void DrawTimeSeriesPlotPlane(void);		//!< 時系列プロット平面を描画する関数
		void DrawTimeSeriesPlot(void);			//!< 時系列プロットを描画する関数
		void DrawWorkSpacePlotPlane(void);		//!< 作業空間プロット平面を描画する関数
		void DrawWorkSpacePlot(void);			//!< 作業空間プロットを描画する関数
		void DrawUserPlotPlane(void);			//!< ユーザプロット平面を描画する関数
		void DrawUserPlot(void);				//!< ユーザプロットを描画する関数
		
		// フレームバッファとキュイプロット
		FrameGraphics FG;						//!< フレームバッファ
		std::array<std::unique_ptr<CuiPlot>, ConstParams::PLOT_NUM> Plot;		//!< 時系列用キュイプロットへのスマートポインタのクラス配列
		CuiPlot PlotXY;							//!< XY作業空間用キュイプロット
		CuiPlot PlotXZ;							//!< XZ作業空間用キュイプロット
		
		// 時系列プロット読み込み用変数
		pthread_mutex_t PlotVarsMutex;	//!< プロット描画変数用のMutex
		bool StorageEnable;				//!< リングバッファ保存時刻になったかの判定用フラグ
		size_t PlotNumBuf;				//!< プロット平面番号バッファ
		size_t VarsCount;				//!< 再帰カウンタ
		
		// 時系列用リングバッファ
		RingBuffer<double, ConstParams::PLOT_RINGBUFF, false> TimeRingBuf;		//!< 時間リングバッファ
		std::array<
			std::array<
				RingBuffer<double, ConstParams::PLOT_RINGBUFF, false>,
				ConstParams::PLOT_VAR_MAX
			>,
			ConstParams::PLOT_NUM
		> VarsRingBuf;		//!< 変数値リングバッファの2次元配列
		
		// 作業空間用バッファ
		pthread_mutex_t WorkspaceMutex;		//!< 作業空間用のMutex
		std::array<Matrix<1,6>, 6> AxisPos;	//!< [m,m,m,0,0,0] XYZ--- 1軸～6軸の作業空間位置
};
}

#endif

