//! @file CustomPlot.hh
//! @brief カスタムプロットクラス
//!
//! ユーザが自由にカスタマイズできるグラフプロット描画クラス
//!
//! @date 2024/06/19
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef CUSTOMPLOT
#define CUSTOMPLOT

#include <cassert>
#include "ARCSassert.hh"
#include "ARCSeventlog.hh"
#include "FrameGraphics.hh"
#include "CuiPlot.hh"

namespace ARCS {	// ARCS名前空間
//! @brief カスタムプロットクラス
class CustomPlot {
	public:
		//! @brief コンストラクタ
		CustomPlot(FrameGraphics& Frame)
			: FG(Frame),	// フレームバッファへの参照を格納
			  Plot(FG, PLOT_LEFT, PLOT_TOP, PLOT_WIDTH, PLOT_HEIGHT)	// プロットの設定（例）
		{
			PassedLog();
		}

		//! @brief デストラクタ
		~CustomPlot(){
			PassedLog();
		}
		
		//! @brief カスタムプロット平面を描画する関数
		void DrawPlotPlane(void){
			// カスタムプロットのグラフパラメータの設定＆描画（例）
			Plot.Visible(PLOT_VISIBLE);	// 可視化設定
			Plot.SetColors(
				PLOT_AXIS_COLOR,		// 軸の色の設定
				PLOT_GRID_COLOR,		// グリッドの色の設定
				PLOT_TEXT_COLOR,		// 文字色の設定
				PLOT_BACK_COLOR,		// 背景色の設定
				PLOT_CURS_COLOR			// カーソルの色の設定
			);
			Plot.SetAxisLabels(PLOT_XLABEL, PLOT_YLABEL);				// 軸ラベルの設定
			Plot.SetRanges(PLOT_XMIN, PLOT_XMAX, PLOT_YMIN, PLOT_YMAX);	// 軸の範囲設定
			Plot.SetGridDivision(PLOT_XGRID, PLOT_YGRID);				// グリッドの分割数の設定
			Plot.DrawAxis();											// 軸の描画
			Plot.StorePlaneInBuffer();									// プロット平面の描画データをバッファに保存しておく
			Plot.Disp();												// プロット平面を画面表示
		}
		
		//! @brief カスタムプロットを描画する関数
		void DrawPlot(void){
			// カスタムプロットの描画動作（例）
			//Plot.LoadPlaneFromBuffer();	// 背景のプロット平面をバッファから読み出す
			
				// ここに時間で変動するプロットを記述する
			
			//Plot.Disp();				// プロット平面＋プロットの描画
		}
		
	private:
		CustomPlot(const CustomPlot&) = delete;					//!< コピーコンストラクタ使用禁止
		CustomPlot(CustomPlot&&) = delete;						//!< ムーブコンストラクタ使用禁止
		const CustomPlot& operator=(const CustomPlot&) = delete;//!< 代入演算子使用禁止
		
		// カスタムプロットの設定（例）
		static constexpr bool PLOT_VISIBLE = false;			//!< プロット可視/不可視設定
		static constexpr int PLOT_LEFT = 100;				//!< [px] 左位置
		static constexpr int PLOT_TOP = 100;				//!< [px] 上位置
		static constexpr int PLOT_WIDTH = 300;				//!< [px] 幅
		static constexpr int PLOT_HEIGHT = 270;				//!< [px] 高さ
		static constexpr char PLOT_XLABEL[] = "X AXIS [-]";	//!< X軸ラベル
		static constexpr char PLOT_YLABEL[] = "Y AXIS [-]";	//!< Y軸ラベル
		static constexpr double PLOT_XMAX =  10;			//!< [mm] X軸最大値
		static constexpr double PLOT_XMIN = -10;			//!< [mm] X軸最小値
		static constexpr double PLOT_YMAX =  20;			//!< [mm] Y軸最大値
		static constexpr double PLOT_YMIN =   0;			//!< [mm] Y軸最小値
		static constexpr unsigned int PLOT_XGRID = 4;		//!< X軸グリッドの分割数
		static constexpr unsigned int PLOT_YGRID = 4;		//!< Y軸グリッドの分割数
		static constexpr FGcolors PLOT_AXIS_COLOR = FGcolors::WHITE;	//!< 軸の色
		static constexpr FGcolors PLOT_GRID_COLOR = FGcolors::GRAY25;	//!< グリッドの色
		static constexpr FGcolors PLOT_BACK_COLOR = FGcolors::BLACK;	//!< 背景色
		static constexpr FGcolors PLOT_TEXT_COLOR = FGcolors::WHITE;	//!< 文字色
		static constexpr FGcolors PLOT_CURS_COLOR = FGcolors::GRAY50;	//!< 時刻カーソルの色
		
		FrameGraphics& FG;	//!< フレームグラフィックスへの参照
		CuiPlot Plot;		//!< キュイプロット（例）
};
}

#endif

