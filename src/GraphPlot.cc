//! @file GraphPlot.cc
//! @brief グラフプロットクラス
//!
//! グラフを描画するクラス
//!
//! @date 2021/09/29
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <memory>
#include <cmath>
#include "GraphPlot.hh"

using namespace ARCS;

//! @brief コンストラクタ
GraphPlot::GraphPlot(void)
	: UserPlot(FG, ConstParams::PLOTUS_LEFT, ConstParams::PLOTUS_TOP, ConstParams::PLOTUS_WIDTH, ConstParams::PLOTUS_HEIGHT),
	  FG(ConstParams::PLOT_FRAMEBUFF),
	  Plot({nullptr}),
	  PlotXY(FG, ConstParams::PLOTXY_LEFT, ConstParams::PLOTXY_TOP, ConstParams::PLOTXY_WIDTH, ConstParams::PLOTXY_HEIGHT),
	  PlotXZ(FG, ConstParams::PLOTXZ_LEFT, ConstParams::PLOTXZ_TOP, ConstParams::PLOTXZ_WIDTH, ConstParams::PLOTXZ_HEIGHT),
	  PlotVarsMutex(PTHREAD_MUTEX_INITIALIZER),
	  StorageEnable(false),
	  PlotNumBuf(0),
	  VarsCount(0),
	  TimeRingBuf(),
	  VarsRingBuf(),
	  WorkspaceMutex(PTHREAD_MUTEX_INITIALIZER),
	  AxisPos()
{
	PassedLog();
	
	// Mutex初期化
	pthread_mutex_init(&PlotVarsMutex, nullptr);
	pthread_mutex_init(&WorkspaceMutex, nullptr);
	
	// 時系列プロット平面の分だけキュイプロットを生成
	for(size_t i = 0; i < ConstParams::PLOT_NUM; ++i){
		Plot.at(i) = std::make_unique<CuiPlot>(
			FG, ConstParams::PLOT_LEFT[i], ConstParams::PLOT_TOP[i], ConstParams::PLOT_WIDTH[i], ConstParams::PLOT_HEIGHT[i]
		);
	}
	
	PassedLog();
}

//! @brief デストラクタ
GraphPlot::~GraphPlot(){
	PassedLog();
}

//! @brief 作業空間プロットに位置ベクトルを設定する関数（6軸配列版）
//! @param[in]	AxPosition	1軸～6軸の作業空間位置ベクトル XYZ--- [m,m,m,0,0,0]^T
void GraphPlot::SetWorkspace(const std::array<Matrix<1,6>, 6>& AxPosition){
	pthread_mutex_lock(&WorkspaceMutex);
	AxisPos = AxPosition;
	pthread_mutex_unlock(&WorkspaceMutex);
}

//! @brief 作業空間プロットに位置ベクトルを設定する関数（3軸個別版）
//! @param[in]	Pos1	1軸の作業空間位置ベクトル XYZ--- [m,m,m]^T
//! @param[in]	Pos2	2軸の作業空間位置ベクトル XYZ--- [m,m,m]^T
//! @param[in]	Pos3	3軸の作業空間位置ベクトル XYZ--- [m,m,m]^T
void GraphPlot::SetWorkspace(const Matrix<1,3>& Pos1, const Matrix<1,3>& Pos2, const Matrix<1,3>& Pos3){
	pthread_mutex_lock(&WorkspaceMutex);
	AxisPos.at(3).LoadShortVector(Pos1);
	AxisPos.at(4).LoadShortVector(Pos2);
	AxisPos.at(5).LoadShortVector(Pos3);
	pthread_mutex_unlock(&WorkspaceMutex);
}

//! @brief プロット平面の描画
void GraphPlot::DrawPlotPlane(void){
	DrawTimeSeriesPlotPlane();	// 時系列プロット平面の描画
	DrawWorkSpacePlotPlane();	// 作業空間プロット平面の描画
	DrawUserPlotPlane();		// ユーザプロット平面の描画
}

//! @brief プロット波形の描画
void GraphPlot::DrawWaves(void){
	DrawTimeSeriesPlot();		// 時系列プロットの描画
	DrawWorkSpacePlot();		// 作業空間プロットの描画
	DrawUserPlot();				// ユーザプロットの描画
}

//! @brief 再開始後にプロットをリセットする関数
void GraphPlot::ResetWaves(void){
	TimeRingBuf.ClearBuffer();		// 時間リングバッファをクリア
	
	// 時系列プロット平面の分だけ回す
	for(size_t j = 0; j < ConstParams::PLOT_NUM; ++j){
		// 変数の分ごとの時系列データのプロット
		for(size_t i = 0; i < ConstParams::PLOT_VAR_NUM[j]; ++i){
			VarsRingBuf.at(j).at(i).ClearBuffer();	// 変数値リングバッファをクリア
		}
	}
}

//! @brief 画面をPNGファイルとして出力する関数
void GraphPlot::SaveScreenImage(void){
	EventLog("Writing PNG File...");
	FG.LoadFrameToScreen();								// フレームバッファから画面バッファに読み込み
	FG.SavePngImageFile(ConstParams::PLOT_PNGFILENAME);	// PNGファイル書き出し
	EventLog("Writing PNG File...Done");
}

//! @brief 時系列プロット平面を描画する関数
void GraphPlot::DrawTimeSeriesPlotPlane(void){
	// 時系列プロット平面の分だけグラフパラメータの設定＆描画
	for(size_t j = 0; j < ConstParams::PLOT_NUM; ++j){
		Plot.at(j)->Visible(ConstParams::PLOT_VISIBLE.at(j));	// 可視化設定
		Plot.at(j)->SetColors(
			ConstParams::PLOT_AXIS_COLOR,	// 軸の色の設定
			ConstParams::PLOT_GRID_COLOR,	// グリッドの色の設定
			ConstParams::PLOT_TEXT_COLOR,	// 文字色の設定
			ConstParams::PLOT_BACK_COLOR,	// 背景色の設定
			ConstParams::PLOT_CURS_COLOR	// 時刻カーソルの色の設定
		);
		Plot.at(j)->SetAxisLabels(ConstParams::PLOT_TLABEL, ConstParams::PLOT_FLABEL.at(j));								// 軸ラベルの設定
		Plot.at(j)->SetRanges(0, ConstParams::PLOT_TIMESPAN, ConstParams::PLOT_FMIN.at(j), ConstParams::PLOT_FMAX.at(j));	// 軸の範囲設定
		Plot.at(j)->SetGridDivision(ConstParams::PLOT_TGRID_NUM, ConstParams::PLOT_FGRID_NUM.at(j));						// グリッドの分割数の設定
		Plot.at(j)->SetGridLabelFormat(ConstParams::PLOT_TFORMAT, ConstParams::PLOT_FFORMAT.at(j));							// グリッドのラベルの書式設定
		Plot.at(j)->DrawAxis();				// 軸の描画
		Plot.at(j)->DrawLegends(ConstParams::PLOT_VAR_NAMES.at(j), ConstParams::PLOT_VAR_COLORS, ConstParams::PLOT_VAR_NUM.at(j));	// 凡例の設定＆描画
		Plot.at(j)->StorePlaneInBuffer();	// プロット平面の描画データをバッファに保存しておく
		Plot.at(j)->Disp();					// プロット平面を画面表示
	}
}

//! @brief 時系列プロットを描画する関数
void GraphPlot::DrawTimeSeriesPlot(void){
	// 時系列プロット平面の分だけ描画
	for(size_t j = 0; j < ConstParams::PLOT_NUM; ++j){
		Plot.at(j)->LoadPlaneFromBuffer();	// 背景のプロット平面をバッファから読み出す
		// 変数の分ごとの時系列データのプロット
		for(size_t i = 0; i < ConstParams::PLOT_VAR_NUM.at(j); ++i){
			pthread_mutex_lock(&PlotVarsMutex);
			Plot.at(j)->TimeSeriesPlot(TimeRingBuf, VarsRingBuf.at(j).at(i), ConstParams::PLOT_TYPE.at(j).at(i), ConstParams::PLOT_VAR_COLORS.at(i));
			pthread_mutex_unlock(&PlotVarsMutex);
		}
		Plot.at(j)->Disp();					// プロット平面＋プロットの描画
	}
}

//! @brief 作業空間プロット平面を描画する関数
void GraphPlot::DrawWorkSpacePlotPlane(void){
	// 作業空間XYプロットのグラフパラメータの設定＆描画
	PlotXY.Visible(ConstParams::PLOTXY_VISIBLE);	// 可視化設定
	PlotXY.SetColors(
		ConstParams::PLOT_AXIS_COLOR,	// 軸の色の設定
		ConstParams::PLOT_GRID_COLOR,	// グリッドの色の設定
		ConstParams::PLOT_TEXT_COLOR,	// 文字色の設定
		ConstParams::PLOT_BACK_COLOR,	// 背景色の設定
		ConstParams::PLOT_CURS_COLOR	// カーソルの色の設定
	);
	PlotXY.SetAxisLabels(ConstParams::PLOTXY_XLABEL, ConstParams::PLOTXY_YLABEL);	// 軸ラベルの設定
	PlotXY.SetRanges(ConstParams::PLOTXY_XMIN, ConstParams::PLOTXY_XMAX, ConstParams::PLOTXY_YMIN, ConstParams::PLOTXY_YMAX);	// 軸の範囲設定
	PlotXY.SetGridDivision(ConstParams::PLOTXY_XGRID, ConstParams::PLOTXY_YGRID);	// グリッドの分割数の設定
	PlotXY.DrawAxis();																// 軸の描画
	PlotXY.Plot(0, 0, CuiPlotTypes::PLOT_CROSS, FGcolors::CYAN);					// ロボットベース原点位置の描画
	PlotXY.StorePlaneInBuffer();													// プロット平面の描画データをバッファに保存しておく
	PlotXY.Disp();																	// プロット平面を画面表示
	
	// 作業空間XZプロットのグラフパラメータの設定＆描画
	PlotXZ.Visible(ConstParams::PLOTXZ_VISIBLE);	// 可視化設定
	PlotXZ.SetColors(
		ConstParams::PLOT_AXIS_COLOR,	// 軸の色の設定
		ConstParams::PLOT_GRID_COLOR,	// グリッドの色の設定
		ConstParams::PLOT_TEXT_COLOR,	// 文字色の設定
		ConstParams::PLOT_BACK_COLOR,	// 背景色の設定
		ConstParams::PLOT_CURS_COLOR	// カーソルの色の設定
	);
	PlotXZ.SetAxisLabels(ConstParams::PLOTXZ_XLABEL, ConstParams::PLOTXZ_ZLABEL);	// 軸ラベルの設定
	PlotXZ.SetRanges(ConstParams::PLOTXZ_XMIN, ConstParams::PLOTXZ_XMAX, ConstParams::PLOTXZ_ZMIN, ConstParams::PLOTXZ_ZMAX);	// 軸の範囲設定
	PlotXZ.SetGridDivision(ConstParams::PLOTXZ_XGRID, ConstParams::PLOTXZ_ZGRID);	// グリッドの分割数の設定
	PlotXZ.DrawAxis();																// 軸の描画
	PlotXZ.Plot(0, 0, CuiPlotTypes::PLOT_CROSS, FGcolors::CYAN);					// ロボットベース原点位置の描画
	PlotXZ.StorePlaneInBuffer();													// プロット平面の描画データをバッファに保存しておく
	PlotXZ.Disp();																	// プロット平面＋プロットを画面表示
}

//! @brief 作業空間プロットを描画する関数
void GraphPlot::DrawWorkSpacePlot(void){
	// 作業空間XYプロットの描画
	PlotXY.LoadPlaneFromBuffer();	// 背景のプロット平面をバッファから読み出す
	pthread_mutex_lock(&WorkspaceMutex);
	PlotXY.Plot(               0,                0, AxisPos.at(0)[1], AxisPos.at(0)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// ベースから1軸までの線
	PlotXY.Plot(AxisPos.at(0)[1], AxisPos.at(0)[2], AxisPos.at(1)[1], AxisPos.at(1)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 1軸から2軸までの線
	PlotXY.Plot(AxisPos.at(1)[1], AxisPos.at(1)[2], AxisPos.at(2)[1], AxisPos.at(2)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 2軸から3軸までの線
	PlotXY.Plot(AxisPos.at(2)[1], AxisPos.at(2)[2], AxisPos.at(3)[1], AxisPos.at(3)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 3軸から4軸までの線
	PlotXY.Plot(AxisPos.at(3)[1], AxisPos.at(3)[2], AxisPos.at(4)[1], AxisPos.at(4)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 4軸から5軸までの線
	PlotXY.Plot(AxisPos.at(4)[1], AxisPos.at(4)[2], AxisPos.at(5)[1], AxisPos.at(5)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 5軸から6軸までの線
	PlotXY.Plot(AxisPos.at(0)[1], AxisPos.at(0)[2], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 1軸の点
	PlotXY.Plot(AxisPos.at(1)[1], AxisPos.at(1)[2], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 2軸の点
	PlotXY.Plot(AxisPos.at(2)[1], AxisPos.at(2)[2], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 3軸の点
	PlotXY.Plot(AxisPos.at(3)[1], AxisPos.at(3)[2], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 4軸の点
	PlotXY.Plot(AxisPos.at(4)[1], AxisPos.at(4)[2], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 5軸の点
	PlotXY.Plot(AxisPos.at(5)[1], AxisPos.at(5)[2], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 6軸の点
	PlotXY.DrawValue(ConstParams::PLOTXY_VAL_XPOS, ConstParams::PLOTXY_VAL_YPOS      , "X = % 7.1f mm", AxisPos.at(5)[1]*1e3);	// X位置数値表示
	PlotXY.DrawValue(ConstParams::PLOTXY_VAL_XPOS, ConstParams::PLOTXY_VAL_YPOS - 0.1, "Y = % 7.1f mm", AxisPos.at(5)[2]*1e3);	// Y位置数値表示
	PlotXY.DrawValue(ConstParams::PLOTXY_VAL_XPOS, ConstParams::PLOTXY_VAL_YPOS - 0.2, "Z = % 7.1f mm", AxisPos.at(5)[3]*1e3);	// Z位置数値表示
	pthread_mutex_unlock(&WorkspaceMutex);
	PlotXY.Disp();	// プロット平面＋プロットの描画
	
	// 作業空間XZプロットの描画
	PlotXZ.LoadPlaneFromBuffer();	// 背景のプロット平面をバッファから読み出す
	pthread_mutex_lock(&WorkspaceMutex);
	PlotXZ.Plot(               0,                0, AxisPos.at(0)[1], AxisPos.at(0)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// ベースから1軸までの線
	PlotXZ.Plot(AxisPos.at(0)[1], AxisPos.at(0)[3], AxisPos.at(1)[1], AxisPos.at(1)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 1軸から2軸までの線
	PlotXZ.Plot(AxisPos.at(1)[1], AxisPos.at(1)[3], AxisPos.at(2)[1], AxisPos.at(2)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 2軸から3軸までの線
	PlotXZ.Plot(AxisPos.at(2)[1], AxisPos.at(2)[3], AxisPos.at(3)[1], AxisPos.at(3)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 3軸から4軸までの線
	PlotXZ.Plot(AxisPos.at(3)[1], AxisPos.at(3)[3], AxisPos.at(4)[1], AxisPos.at(4)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 4軸から5軸までの線
	PlotXZ.Plot(AxisPos.at(4)[1], AxisPos.at(4)[3], AxisPos.at(5)[1], AxisPos.at(5)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// 5軸から6軸までの線
	PlotXZ.Plot(AxisPos.at(0)[1], AxisPos.at(0)[3], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 1軸の点
	PlotXZ.Plot(AxisPos.at(1)[1], AxisPos.at(1)[3], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 2軸の点
	PlotXZ.Plot(AxisPos.at(2)[1], AxisPos.at(2)[3], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 3軸の点
	PlotXZ.Plot(AxisPos.at(3)[1], AxisPos.at(3)[3], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 4軸の点
	PlotXZ.Plot(AxisPos.at(4)[1], AxisPos.at(4)[3], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 5軸の点
	PlotXZ.Plot(AxisPos.at(5)[1], AxisPos.at(5)[3], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// 6軸の点
	PlotXZ.DrawValue(ConstParams::PLOTXZ_VAL_XPOS, ConstParams::PLOTXZ_VAL_ZPOS      , "R = % 6.1f deg", AxisPos.at(5)[4]*180.0/M_PI);	// ロール角数値表示
	PlotXZ.DrawValue(ConstParams::PLOTXZ_VAL_XPOS, ConstParams::PLOTXZ_VAL_ZPOS - 0.1, "P = % 6.1f deg", AxisPos.at(5)[5]*180.0/M_PI);	// ピッチ角数値表示
	PlotXZ.DrawValue(ConstParams::PLOTXZ_VAL_XPOS, ConstParams::PLOTXZ_VAL_ZPOS - 0.2, "W = % 6.1f deg", AxisPos.at(5)[6]*180.0/M_PI);	// ヨー角数値表示
	pthread_mutex_unlock(&WorkspaceMutex);
	PlotXZ.Disp();	// プロット平面＋プロットの描画
}

//! @brief ユーザプロット平面を描画する関数
void GraphPlot::DrawUserPlotPlane(void){
	// ユーザプロットのグラフパラメータの設定＆描画
	UserPlot.Visible(ConstParams::PLOTUS_VISIBLE);	// 可視化設定
	UserPlot.SetColors(
		ConstParams::PLOT_AXIS_COLOR,	// 軸の色の設定
		ConstParams::PLOT_GRID_COLOR,	// グリッドの色の設定
		ConstParams::PLOT_TEXT_COLOR,	// 文字色の設定
		ConstParams::PLOT_BACK_COLOR,	// 背景色の設定
		ConstParams::PLOT_CURS_COLOR	// カーソルの色の設定
	);
	UserPlot.SetAxisLabels(ConstParams::PLOTUS_XLABEL, ConstParams::PLOTUS_YLABEL);	// 軸ラベルの設定
	UserPlot.SetRanges(ConstParams::PLOTUS_XMIN, ConstParams::PLOTUS_XMAX, ConstParams::PLOTUS_YMIN, ConstParams::PLOTUS_YMAX);	// 軸の範囲設定
	UserPlot.SetGridDivision(ConstParams::PLOTUS_XGRID, ConstParams::PLOTUS_YGRID);	// グリッドの分割数の設定
	UserPlot.DrawAxis();															// 軸の描画
	UserPlot.StorePlaneInBuffer();													// プロット平面の描画データをバッファに保存しておく
	UserPlot.Disp();																// プロット平面を画面表示
}

//! @brief ユーザプロットを描画する関数
void GraphPlot::DrawUserPlot(void){
	//PlotUser.LoadPlaneFromBuffer();	// 背景のプロット平面をバッファから読み出す
	
		// ここに時間で変動するプロットを記述する
	
	//PlotUser.Disp();				// プロット平面＋プロットの描画
}

