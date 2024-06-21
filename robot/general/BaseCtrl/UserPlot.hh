//! @file UserPlot.hh
//! @brief ユーザカスタムプロットクラス
//!
//! ユーザが自由にカスタマイズできるグラフプロット描画クラス
//!
//! @date 2024/06/21
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef USERPLOT
#define USERPLOT

#include <cassert>
#include <functional>
#include "ARCSassert.hh"
#include "ARCSeventlog.hh"
#include "GraphPlot.hh"
#include "CuiPlot.hh"

namespace ARCS {	// ARCS名前空間
//! @brief ユーザカスタムプロットクラス
class UserPlot {
	public:
		//! @brief コンストラクタ
		UserPlot(GraphPlot& GP)
			: // 以下でグラフ描画に使用するクラスと変数を初期化
			  Plot(GP.GetFGrefs() , PLOT_LEFT, PLOT_TOP, PLOT_WIDTH, PLOT_HEIGHT),	// キュイプロットの設定（例）
			  X1(0), Y1(0),	// プロット変数（例）

			  // 以下は編集しないこと
			  Graph(GP), DrawPlaneFobj(), DrawPlotFobj()	// 初期化子
		{
			Initialize();	// 初期化
		}
		
		//! @brief プロット変数設定関数（例）
		void SetVars(const double x, const double y){
			// 以下にグラフ描画したい変数値を設定（例）
			X1 = x;
			Y1 = y;
		}

	private:
		// ユーザカスタムプロットの設定（例）
		static constexpr bool PLOT_VISIBLE = true;			//!< プロット可視/不可視設定
		static constexpr int PLOT_LEFT = 1015;				//!< [px] 左位置
		static constexpr int PLOT_TOP = 97;					//!< [px] 上位置
		static constexpr int PLOT_WIDTH = 300;				//!< [px] 幅
		static constexpr int PLOT_HEIGHT = 270;				//!< [px] 高さ
		static constexpr char PLOT_XLABEL[] = "X AXIS [-]";	//!< X軸ラベル
		static constexpr char PLOT_YLABEL[] = "Y AXIS [-]";	//!< Y軸ラベル
		static constexpr double PLOT_XMAX =  10;			//!< [mm] X軸最大値
		static constexpr double PLOT_XMIN = -10;			//!< [mm] X軸最小値
		static constexpr double PLOT_YMAX =  10;			//!< [mm] Y軸最大値
		static constexpr double PLOT_YMIN = -10;			//!< [mm] Y軸最小値
		static constexpr unsigned int PLOT_XGRID = 4;		//!< X軸グリッドの分割数
		static constexpr unsigned int PLOT_YGRID = 4;		//!< Y軸グリッドの分割数
		static constexpr FGcolors PLOT_AXIS_COLOR = FGcolors::WHITE;	//!< 軸の色
		static constexpr FGcolors PLOT_GRID_COLOR = FGcolors::GRAY25;	//!< グリッドの色
		static constexpr FGcolors PLOT_BACK_COLOR = FGcolors::BLACK;	//!< 背景色
		static constexpr FGcolors PLOT_TEXT_COLOR = FGcolors::WHITE;	//!< 文字色
		static constexpr FGcolors PLOT_CURS_COLOR = FGcolors::GRAY50;	//!< 時刻カーソルの色
		
		// 以下にグラフ描画に使用するクラスと変数を定義
		CuiPlot Plot;		//!< キュイプロット（例）
		double X1, Y1;		//!< プロット変数（例）
		
		//! @brief ユーザカスタムプロット平面を描画する関数
		//! （この関数はGraphPlotクラスの内部で関数オブジェクトを介して開始時に一度だけ実行される）
		void DrawPlotPlane(void){
			// ユーザカスタムプロットのグラフパラメータの設定＆描画（例）
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
		
		//! @brief ユーザカスタムプロットを描画する関数
		//! （この関数はGraphPlotクラスの内部で関数オブジェクトを介して描画周期毎に実行される）
		void DrawPlot(void){
			// ユーザカスタムプロットの描画動作（例）
			//Plot.LoadPlaneFromBuffer();	// 背景のプロット平面をバッファから読み出す
			
			// ここに時間で変動するプロットを記述する
			Plot.Plot(X1, Y1, CuiPlotTypes::PLOT_CROSS, FGcolors::ORANGE);	// データ点を1点プロット（例）
			
			Plot.Disp();	// プロット平面＋プロットの描画
		}

	// ここから下は編集しないこと
	public:
		//! @brief デストラクタ
		~UserPlot(){
			PassedLog();
		}

	private:
		UserPlot(const UserPlot&) = delete;					//!< コピーコンストラクタ使用禁止
		UserPlot(UserPlot&&) = delete;						//!< ムーブコンストラクタ使用禁止
		const UserPlot& operator=(const UserPlot&) = delete;//!< 代入演算子使用禁止
		GraphPlot& Graph;									//!< グラフプロットへの参照
		std::function<void(void)> DrawPlaneFobj;			//!< ユーザカスタムプロット平面描画関数の関数オブジェクト
		std::function<void(void)> DrawPlotFobj;				//!< ユーザカスタムプロット描画関数の関数オブジェクト

		//! @brief 初期化関数
		void Initialize(void){
			PassedLog();
			DrawPlaneFobj = [&](void){ return DrawPlotPlane(); };	// プロット平面描画関数への関数オブジェクトをラムダ式で格納
			DrawPlotFobj  = [&](void){ return DrawPlot(); };		// プロット描画関数への関数オブジェクトをラムダ式で格納
			Graph.SetUserPlotFuncs(DrawPlaneFobj, DrawPlotFobj);	// 描画関数オブジェクトをグラフプロットへ渡す
		}
};
}

#endif

