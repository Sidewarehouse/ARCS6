//! @file ARCSgraphics.cc
//! @brief グラフィッククラス
//!
//! グラフを描画するクラス
//!
//! @date 2024/10/08
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <memory>
#include <cmath>
#include "ARCSgraphics.hh"
#include "EquipParams.hh"

using namespace ARCS;

//! @brief コンストラクタ
ARCSgraphics::ARCSgraphics(void)
	: FG(EquipParams::PLOT_FRAMEBUFF),
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
	  AxisPos(),
	  DrawUserPlaneFunc(),	
	  DrawUserPlotFunc()
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
ARCSgraphics::~ARCSgraphics(){
	PassedLog();
}

//! @brief 作業空間プロットに位置ベクトルを設定する関数
//! @param[in]	AxPosition	1軸～6軸の作業空間位置ベクトル XYZ--- [m,m,m,0,0,0]^T
void ARCSgraphics::SetWorkspace(const std::array<ArcsMat<6,1>, ConstParams::PLOTXYXZ_NUMPT>& AxPosition){
	//pthread_mutex_lock(&WorkspaceMutex);	 ←リアルタイム性悪化の原因！(一時的な対処)
	AxisPos = AxPosition;
	//pthread_mutex_unlock(&WorkspaceMutex);
}

//! @brief フレームバッファクラスへの参照を返す関数
//! @return	FrameGraphicsへの参照
FrameGraphics<>& ARCSgraphics::GetFGrefs(void){
	return FG;
}

//! @brief ユーザカスタムプロット描画関数への関数オブジェクトを設定する関数
void ARCSgraphics::SetUserPlotFuncs(const std::function<void(void)>& DrawPlaneFobj, const std::function<void(void)>& DrawPlotFobj){
	DrawUserPlaneFunc = DrawPlaneFobj;	// ユーザカスタムプロット平面描画関数の関数オブジェクトを格納
	DrawUserPlotFunc  = DrawPlotFobj;	// ユーザカスタムプロット描画関数の関数オブジェクトを格納
}

//! @brief プロット平面の描画
void ARCSgraphics::DrawPlotPlane(void){
	DrawTimeSeriesPlotPlane();	// 時系列プロット平面の描画
	DrawWorkSpacePlotPlane();	// 作業空間プロット平面の描画

	// 関数オブジェクトが準備完了するまで待機（同期機構を入れるべきだが暫定措置）
	while(!DrawUserPlaneFunc){
		usleep(ARCSparams::ARCS_TIME_GRPL);
	}
	DrawUserPlaneFunc(); 		// ユーザカスタムプロット平面の描画(関数オブジェクト経由で実行)
}

//! @brief プロット波形の描画
void ARCSgraphics::DrawWaves(void){
	DrawTimeSeriesPlot();		// 時系列プロットの描画
	DrawWorkSpacePlot();		// 作業空間プロットの描画
	DrawUserPlotFunc();			// ユーザカスタムプロットの描画(関数オブジェクト経由で実行)
}

//! @brief 再開始後にプロットをリセットする関数
void ARCSgraphics::ResetWaves(void){
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
void ARCSgraphics::SaveScreenImage(void){
	EventLog("Writing PNG File...");
	FG.LoadFrameToScreen();								// フレームバッファから画面バッファに読み込み
	FG.SavePngImageFile(ConstParams::PLOT_PNGFILENAME);	// PNGファイル書き出し
	EventLog("Writing PNG File...Done");
}

//! @brief 時系列プロット平面を描画する関数
void ARCSgraphics::DrawTimeSeriesPlotPlane(void){
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
		
		// プロット変数名の char[] を std::string に変換
		std::array<std::string, ARCSparams::PLOT_VAR_MAX> VarNameBuff;
		for(size_t i = 0; i < ARCSparams::PLOT_VAR_MAX; ++i) VarNameBuff.at(i) = ConstParams::PLOT_VAR_NAMES.at(j).at(i);
		
		Plot.at(j)->DrawLegends(VarNameBuff, ConstParams::PLOT_VAR_COLORS, ConstParams::PLOT_VAR_NUM.at(j));	// 凡例の設定＆描画
		Plot.at(j)->StorePlaneInBuffer();	// プロット平面の描画データをバッファに保存しておく
		Plot.at(j)->Disp();					// プロット平面を画面表示
	}
}

//! @brief 時系列プロットを描画する関数
void ARCSgraphics::DrawTimeSeriesPlot(void){
	// 時系列プロット平面の分だけ描画
	for(size_t j = 0; j < ConstParams::PLOT_NUM; ++j){
		Plot.at(j)->LoadPlaneFromBuffer();	// 背景のプロット平面をバッファから読み出す
		// 変数の分ごとの時系列データのプロット
		for(size_t i = 0; i < ConstParams::PLOT_VAR_NUM.at(j); ++i){
			//pthread_mutex_lock(&PlotVarsMutex); ←リアルタイム性悪化の原因！(一時的な対処)
			Plot.at(j)->TimeSeriesPlot(TimeRingBuf, VarsRingBuf.at(j).at(i), ConstParams::PLOT_TYPE.at(j).at(i), ConstParams::PLOT_VAR_COLORS.at(i));
			//pthread_mutex_unlock(&PlotVarsMutex);
		}
		Plot.at(j)->Disp();					// プロット平面＋プロットの描画
	}
}

//! @brief 作業空間プロット平面を描画する関数
void ARCSgraphics::DrawWorkSpacePlotPlane(void){
	// 作業空間XYプロットのグラフパラメータの設定＆描画
	PlotXY.Visible(ConstParams::PLOTXYXZ_VISIBLE);	// 可視化設定
	PlotXY.SetColors(
		ConstParams::PLOT_AXIS_COLOR,	// 軸の色の設定
		ConstParams::PLOT_GRID_COLOR,	// グリッドの色の設定
		ConstParams::PLOT_TEXT_COLOR,	// 文字色の設定
		ConstParams::PLOT_BACK_COLOR,	// 背景色の設定
		ConstParams::PLOT_CURS_COLOR	// カーソルの色の設定
	);
	PlotXY.SetAxisLabels(ConstParams::PLOTXYXZ_XLABEL, ConstParams::PLOTXY_YLABEL);	// 軸ラベルの設定
	PlotXY.SetRanges(ConstParams::PLOTXYXZ_XMIN, ConstParams::PLOTXYXZ_XMAX, ConstParams::PLOTXY_YMIN, ConstParams::PLOTXY_YMAX);	// 軸の範囲設定
	PlotXY.SetGridDivision(ConstParams::PLOTXYXZ_XGRID, ConstParams::PLOTXY_YGRID);	// グリッドの分割数の設定
	PlotXY.DrawAxis();																// 軸の描画
	PlotXY.Plot(0, 0, CuiPlotTypes::PLOT_CROSS, FGcolors::CYAN);					// ロボットベース原点位置の描画
	PlotXY.StorePlaneInBuffer();													// プロット平面の描画データをバッファに保存しておく
	PlotXY.Disp();																	// プロット平面を画面表示
	
	// 作業空間XZプロットのグラフパラメータの設定＆描画
	PlotXZ.Visible(ConstParams::PLOTXYXZ_VISIBLE);	// 可視化設定
	PlotXZ.SetColors(
		ConstParams::PLOT_AXIS_COLOR,	// 軸の色の設定
		ConstParams::PLOT_GRID_COLOR,	// グリッドの色の設定
		ConstParams::PLOT_TEXT_COLOR,	// 文字色の設定
		ConstParams::PLOT_BACK_COLOR,	// 背景色の設定
		ConstParams::PLOT_CURS_COLOR	// カーソルの色の設定
	);
	PlotXZ.SetAxisLabels(ConstParams::PLOTXYXZ_XLABEL, ConstParams::PLOTXZ_ZLABEL);	// 軸ラベルの設定
	PlotXZ.SetRanges(ConstParams::PLOTXYXZ_XMIN, ConstParams::PLOTXYXZ_XMAX, ConstParams::PLOTXZ_ZMIN, ConstParams::PLOTXZ_ZMAX);	// 軸の範囲設定
	PlotXZ.SetGridDivision(ConstParams::PLOTXYXZ_XGRID, ConstParams::PLOTXZ_ZGRID);	// グリッドの分割数の設定
	PlotXZ.DrawAxis();																// 軸の描画
	PlotXZ.Plot(0, 0, CuiPlotTypes::PLOT_CROSS, FGcolors::CYAN);					// ロボットベース原点位置の描画
	PlotXZ.StorePlaneInBuffer();													// プロット平面の描画データをバッファに保存しておく
	PlotXZ.Disp();																	// プロット平面＋プロットを画面表示
}

//! @brief 作業空間プロットを描画する関数
void ARCSgraphics::DrawWorkSpacePlot(void){
	// 作業空間プロット平面の背景をバッファから読み出す
	PlotXY.LoadPlaneFromBuffer();
	PlotXZ.LoadPlaneFromBuffer();

	// 作業空間XYプロットにおける各軸要素の間の線の描画
	// [1] = X軸, [2] = Y軸
	PlotXY.Plot(0, 0, AxisPos.at(0)[1], AxisPos.at(0)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// ベース原点から1軸要素までの線
	for(size_t i = 1; i < ConstParams::PLOTXYXZ_NUMPT; ++i){
		PlotXY.Plot(AxisPos.at(i - 1)[1], AxisPos.at(i - 1)[2], AxisPos.at(i)[1], AxisPos.at(i)[2], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// i-1軸要素からi軸要素までの線
	}

	// 作業空間XYプロットにおける各軸要素の点の描画
	for(size_t i = 0; i < ConstParams::PLOTXYXZ_NUMPT; ++i){
		PlotXY.Plot(AxisPos.at(i)[1], AxisPos.at(i)[2], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);	// i軸要素の点
	}

	// 作業空間XYプロットにおける最先端軸要素の数値表示
	PlotXY.DrawValue(ConstParams::PLOTXY_VAL_XPOS, ConstParams::PLOTXY_VAL_YPOS      , "X = % 7.1f mm", AxisPos.at( ConstParams::PLOTXYXZ_NUMPT - 1 )[1]*1e3);	// X位置数値表示
	PlotXY.DrawValue(ConstParams::PLOTXY_VAL_XPOS, ConstParams::PLOTXY_VAL_YPOS - 0.1, "Y = % 7.1f mm", AxisPos.at( ConstParams::PLOTXYXZ_NUMPT - 1 )[2]*1e3);	// Y位置数値表示
	PlotXY.DrawValue(ConstParams::PLOTXY_VAL_XPOS, ConstParams::PLOTXY_VAL_YPOS - 0.2, "Z = % 7.1f mm", AxisPos.at( ConstParams::PLOTXYXZ_NUMPT - 1 )[3]*1e3);	// Z位置数値表示

	// 作業空間XZプロットにおける各軸要素の間の線の描画
	// [1] = X軸, [3] = Z軸
	PlotXZ.Plot(0, 0, AxisPos.at(0)[1], AxisPos.at(0)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// ベース原点から1軸要素までの線
	for(size_t i = 1; i < ConstParams::PLOTXYXZ_NUMPT; ++i){
		PlotXZ.Plot(AxisPos.at(i - 1)[1], AxisPos.at(i - 1)[3], AxisPos.at(i)[1], AxisPos.at(i)[3], CuiPlotTypes::PLOT_BOLDLINE, FGcolors::CYAN);	// i-1軸要素からi軸要素までの線
	}

	// 作業空間XZプロットにおける各軸要素の点の描画
	for(size_t i = 0; i < ConstParams::PLOTXYXZ_NUMPT; ++i){
		PlotXZ.Plot(AxisPos.at(i)[1], AxisPos.at(i)[3], CuiPlotTypes::PLOT_BOLDDOT, FGcolors::YELLOW);		// i軸要素の点
	}

	// 作業空間XZプロットにおける最先端軸要素の数値表示
	PlotXZ.DrawValue(ConstParams::PLOTXZ_VAL_XPOS, ConstParams::PLOTXZ_VAL_ZPOS      , "R = % 6.1f deg", AxisPos.at( ConstParams::PLOTXYXZ_NUMPT - 1 )[4]*180.0/M_PI);	// ロール角数値表示
	PlotXZ.DrawValue(ConstParams::PLOTXZ_VAL_XPOS, ConstParams::PLOTXZ_VAL_ZPOS - 0.1, "P = % 6.1f deg", AxisPos.at( ConstParams::PLOTXYXZ_NUMPT - 1 )[5]*180.0/M_PI);	// ピッチ角数値表示
	PlotXZ.DrawValue(ConstParams::PLOTXZ_VAL_XPOS, ConstParams::PLOTXZ_VAL_ZPOS - 0.2, "W = % 6.1f deg", AxisPos.at( ConstParams::PLOTXYXZ_NUMPT - 1 )[6]*180.0/M_PI);	// ヨー角数値表示

	// 作業空間プロット平面＋プロットの描画実行
	PlotXY.Disp();
	PlotXZ.Disp();
}

