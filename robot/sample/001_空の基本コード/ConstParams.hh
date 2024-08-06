//! @file ConstParams.hh
//! @brief 定数値格納用クラス
//!        ARCSに必要な定数値を格納します。
//! @date 2024/06/24
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef CONSTPARAMS
#define CONSTPARAMS

#include <cmath>
#include "ARCSparams.hh"
#include "SFthread.hh"
#include "FrameGraphics.hh"
#include "CuiPlot.hh"

namespace ARCS {	// ARCS名前空間
//! @brief 定数値格納用クラス
class ConstParams {
	public:
		// タイトルに表示させる制御系の名前(識別用に好きな名前を入力)
		static constexpr char CTRLNAME[] = "<TITLE: NOTITLE >";	//!< (画面に入る文字数以内)
		
		// 実験データCSVファイルの設定
		static constexpr char DATA_NAME[] = "DATA.csv";	//!< CSVファイル名
		static constexpr double DATA_START =  0;		//!< [s] 保存開始時刻
		static constexpr double DATA_END   = 10;		//!< [s] 保存終了時刻
		static constexpr double DATA_RESO  = 0.001;		//!< [s] データの時間分解能
		static constexpr size_t DATA_NUM  =  10;		//!< [-] 保存する変数の数

		// SCHED_FIFOリアルタイムスレッドの設定
		static constexpr size_t THREAD_NUM = 1;			//!< 動作させるスレッドの数 (最大数は ARCSparams::THREAD_NUM_MAX 個まで)
		
		//! @brief 制御周期の設定
		static constexpr std::array<unsigned long, ARCSparams::THREAD_MAX> SAMPLING_TIME = {
		//   s  m  u  n	制御周期は Ts[0] ≦ Ts[1] ≦ … ≦ Ts[THREAD_MAX] になるようにすること
				 100000,	// [ns] 制御用周期実行関数1 (スレッド1) 制御周期
				1000000,	// [ns] 制御用周期実行関数2 (スレッド2) 制御周期
				1000000,	// [ns] 制御用周期実行関数3 (スレッド3) 制御周期
		};
		
		// デバッグプリントとデバッグインジケータの設定
		static constexpr bool DEBUG_PRINT_VISIBLE = false;	//!< デバッグプリント表示の有効/無効設定
		static constexpr bool DEBUG_INDIC_VISIBLE = false;	//!< デバッグインジケータ表示の有効/無効設定
		
		// 任意変数値表示の設定
		static constexpr size_t INDICVARS_NUM = 10;			//!< 表示したい変数の数 (最大数 INDICVARS_MAX まで)
		
		//! @brief 任意に表示したい変数値の表示形式 (printfの書式と同一)
		static constexpr std::array<char[15], ARCSparams::INDICVARS_MAX> INDICVARS_FORMS = {
			"% 13.4f",	// 変数 0
			"% 13.4f",	// 変数 1
			"% 13.4f",	// 変数 2
			"% 13.4f",	// 変数 3
			"% 13.4f",	// 変数 4
			"% 13.4f",	// 変数 5
			"% 13.4f",	// 変数 6
			"% 13.4f",	// 変数 7
			"% 13.4f",	// 変数 8
			"% 13.4f",	// 変数 9
			"% 13.4f",	// 変数10
			"% 13.4f",	// 変数11
			"% 13.4f",	// 変数12
			"% 13.4f",	// 変数13
			"% 13.4f",	// 変数14
			"% 13.4f",	// 変数15
		};
		
		// オンライン設定変数の設定
		static constexpr size_t ONLINEVARS_NUM = 10;	//!< オンライン設定変数の数 (最大数 ONLINEVARS_MAX まで)
		
		// 時系列グラフプロットの共通設定
		static constexpr char PLOT_PNGFILENAME[] = "Screenshot.png";//!< スクリーンショットのPNGファイル名
		static constexpr size_t PLOT_NUM =  4;						//!< [-] グラフプロットの数
		static constexpr double PLOT_TIMESPAN = 10;					//!< [s] プロットの時間幅
		static constexpr double PLOT_TIMERESO = 0.01;				//!< [s] プロットの時間分解能
		static constexpr size_t PLOT_RINGBUFF = 1024;				//!< [-] プロット用リングバッファの要素数
		static constexpr size_t PLOT_TGRID_NUM = 10;				//!< [-] 時間軸グリッドの分割数
		static constexpr char PLOT_TFORMAT[] = "%3.1f";				//!< 時間軸書式
		static constexpr char PLOT_TLABEL[] = "Time [s]";			//!< 時間軸ラベル
		
		//! @brief 縦軸ラベル
		static constexpr std::array<char[31], ARCSparams::PLOT_MAX> PLOT_FLABEL = {
			"---------- [-]",	// グラフプロット0
			"---------- [-]",	// グラフプロット1
			"---------- [-]",	// グラフプロット2
			"---------- [-]",	// グラフプロット3
			"---------- [-]",	// グラフプロット4
			"---------- [-]",	// グラフプロット5
			"---------- [-]",	// グラフプロット6
			"---------- [-]",	// グラフプロット7
			"---------- [-]",	// グラフプロット8
			"---------- [-]",	// グラフプロット9
			"---------- [-]",	// グラフプロット10
			"---------- [-]",	// グラフプロット11
			"---------- [-]",	// グラフプロット12
			"---------- [-]",	// グラフプロット13
			"---------- [-]",	// グラフプロット14
			"---------- [-]",	// グラフプロット15
		};
		
		//! @brief 縦軸書式
		static constexpr std::array<char[15], ARCSparams::PLOT_MAX> PLOT_FFORMAT = {
			"%6.1f",	// グラフプロット0
			"%6.1f",	// グラフプロット1
			"%6.1f",	// グラフプロット2
			"%6.1f",	// グラフプロット3
			"%6.1f",	// グラフプロット4
			"%6.1f",	// グラフプロット5
			"%6.1f",	// グラフプロット6
			"%6.1f",	// グラフプロット7
			"%6.1f",	// グラフプロット8
			"%6.1f",	// グラフプロット9
			"%6.1f",	// グラフプロット10
			"%6.1f",	// グラフプロット11
			"%6.1f",	// グラフプロット12
			"%6.1f",	// グラフプロット13
			"%6.1f",	// グラフプロット14
			"%6.1f",	// グラフプロット15
		};
		
		//! @brief プロット変数の名前
		static constexpr std::array<
			std::array<char[15], ARCSparams::PLOT_VAR_MAX>, ARCSparams::PLOT_MAX
		> PLOT_VAR_NAMES = {{
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット0
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット1
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット2
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット3
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット4
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット5
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット6
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット7
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット8
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット9
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット10
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット11
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット12
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット13
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット14
			{"VAR-00", "VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット15
		}};
		
		static constexpr FGcolors PLOT_AXIS_COLOR = FGcolors::WHITE;	//!< 軸の色
		static constexpr FGcolors PLOT_GRID_COLOR = FGcolors::GRAY25;	//!< グリッドの色
		static constexpr FGcolors PLOT_BACK_COLOR = FGcolors::BLACK;	//!< 背景色
		static constexpr FGcolors PLOT_TEXT_COLOR = FGcolors::WHITE;	//!< 文字色
		static constexpr FGcolors PLOT_CURS_COLOR = FGcolors::GRAY50;	//!< 時刻カーソルの色
		
		//! @brief 時系列グラフ描画の有効/無効設定
		static constexpr std::array<bool, ARCSparams::PLOT_MAX> PLOT_VISIBLE = {
			true,	// プロット0
			true,	// プロット1
			true,	// プロット2
			true,	// プロット3
			true,	// プロット4
			true,	// プロット5
			true,	// プロット6
			true,	// プロット7
			true,	// プロット8
			true,	// プロット9
			true,	// プロット10
			true,	// プロット11
			true,	// プロット12
			true,	// プロット13
			true,	// プロット14
			true,	// プロット15
		};
		
		//! @brief 時系列プロットの変数ごとの線の色
		static constexpr std::array<FGcolors, ARCSparams::PLOT_VAR_MAX> PLOT_VAR_COLORS = {
			FGcolors::RED,
			FGcolors::GREEN,
			FGcolors::CYAN,
			FGcolors::MAGENTA,
			FGcolors::YELLOW,
			FGcolors::ORANGE,
			FGcolors::WHITE,
			FGcolors::BLUE,
		};
		
		//! @brief 時系列プロットする変数の数 (≦PLOT_VAR_MAX)
		static constexpr std::array<size_t, ARCSparams::PLOT_MAX> PLOT_VAR_NUM = {
			1,	// プロット0
			1,	// プロット1
			1,	// プロット2
			1,	// プロット3
			1,	// プロット4
			1,	// プロット5
			1,	// プロット6
			1,	// プロット7
			1,	// プロット8
			1,	// プロット9
			1,	// プロット10
			1,	// プロット11
			1,	// プロット12
			1,	// プロット13
			1,	// プロット14
			1,	// プロット15
		};
		
		//! @brief 時系列プロットの縦軸最大値
		static constexpr std::array<double, ARCSparams::PLOT_MAX> PLOT_FMAX	= {
			1.0,	// プロット0
			1.0,	// プロット1
			1.0,	// プロット2
			1.0,	// プロット3
			1.0,	// プロット4
			1.0,	// プロット5
			1.0,	// プロット6
			1.0,	// プロット7
			1.0,	// プロット8
			1.0,	// プロット9
			1.0,	// プロット10
			1.0,	// プロット11
			1.0,	// プロット12
			1.0,	// プロット13
			1.0,	// プロット14
			1.0,	// プロット15
		};
		
		//! @brief 時系列プロットの縦軸最小値
		static constexpr std::array<double, ARCSparams::PLOT_MAX> PLOT_FMIN = {
			-1.0,	// プロット0
			-1.0,	// プロット1
			-1.0,	// プロット2
			-1.0,	// プロット3
			-1.0,	// プロット4
			-1.0,	// プロット5
			-1.0,	// プロット6
			-1.0,	// プロット7
			-1.0,	// プロット8
			-1.0,	// プロット9
			-1.0,	// プロット10
			-1.0,	// プロット11
			-1.0,	// プロット12
			-1.0,	// プロット13
			-1.0,	// プロット14
			-1.0,	// プロット15
		};
		
		//! @brief 時系列プロットの縦軸グリッドの分割数
		static constexpr std::array<size_t, ARCSparams::PLOT_MAX> PLOT_FGRID_NUM = {
			4,	// プロット0
			4,	// プロット1
			4,	// プロット2
			4,	// プロット3
			4,	// プロット4
			4,	// プロット5
			4,	// プロット6
			4,	// プロット7
			4,	// プロット8
			4,	// プロット9
			4,	// プロット10
			4,	// プロット11
			4,	// プロット12
			4,	// プロット13
			4,	// プロット14
			4,	// プロット15
		};
		
		//! @brief [px] 時系列プロットの左位置
		static constexpr std::array<int, ARCSparams::PLOT_MAX> PLOT_LEFT = {
			305,	// プロット0
			305,	// プロット1
			305,	// プロット2
			305,	// プロット3
			305,	// プロット4
			305,	// プロット5
			1015,	// プロット6
			1015,	// プロット7
			1015,	// プロット8
			1015,	// プロット9
			1015,	// プロット10
			1015,	// プロット11
			   0,	// プロット12
			   0,	// プロット13
			   0,	// プロット14
			   0,	// プロット15
		};
		
		//! @brief [px] 時系列プロットの上位置
		static constexpr std::array<int, ARCSparams::PLOT_MAX> PLOT_TOP = {
			 97,	// プロット0
			250,	// プロット1
			403,	// プロット2
			556,	// プロット3
			709,	// プロット4
			862,	// プロット5
			 97,	// プロット6
			250,	// プロット7
			403,	// プロット8
			556,	// プロット9
			709,	// プロット10
			862,	// プロット11
			  0,	// プロット12
			  0,	// プロット13
			  0,	// プロット14
			  0,	// プロット15
		};
		
		//! @brief [px] 時系列プロットの幅
		static constexpr std::array<int, ARCSparams::PLOT_MAX> PLOT_WIDTH = {
			710,	// プロット0
			710,	// プロット1
			710,	// プロット2
			710,	// プロット3
			710,	// プロット4
			710,	// プロット5
			710,	// プロット6
			710,	// プロット7
			710,	// プロット8
			710,	// プロット9
			710,	// プロット10
			710,	// プロット11
			710,	// プロット12
			710,	// プロット13
			710,	// プロット14
			710,	// プロット15
		};
		
		//! @brief [px] 時系列プロットの高さ
		static constexpr std::array<int, ARCSparams::PLOT_MAX> PLOT_HEIGHT = {
			153,	// プロット0
			153,	// プロット1
			153,	// プロット2
			153,	// プロット3
			153,	// プロット4
			153,	// プロット5
			153,	// プロット6
			153,	// プロット7
			153,	// プロット8
			153,	// プロット9
			153,	// プロット10
			153,	// プロット11
			153,	// プロット12
			153,	// プロット13
			153,	// プロット14
			153,	// プロット15
		};
		
		//! @brief 時系列プロットの種類の設定
		//! 下記のプロット方法が使用可能
		//!	PLOT_LINE		線プロット
		//!	PLOT_BOLDLINE 	太線プロット
		//!	PLOT_DOT		点プロット
		//!	PLOT_BOLDDOT	太点プロット
		//!	PLOT_CROSS		十字プロット
		//!	PLOT_STAIRS		階段プロット
		//!	PLOT_BOLDSTAIRS	太線階段プロット
		//!	PLOT_LINEANDDOT	線と点の複合プロット
		static constexpr std::array<
			std::array<CuiPlotTypes, ARCSparams::PLOT_VAR_MAX>, ARCSparams::PLOT_MAX
		> PLOT_TYPE = {{
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット0
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット1
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット2
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット3
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット4
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット5
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット6
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット7
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット8
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット9
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット10
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット11
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット12
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット13
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット14
				
			{CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,
				CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE, CuiPlotTypes::PLOT_LINE,},	// プロット15
		}};
		
		//! @brief 作業空間プロット共通の設定
		static constexpr bool PLOTXYXZ_VISIBLE = false;	//!< プロット可視/不可視設定
		static constexpr size_t PLOTXYXZ_NUMPT = 9;		//!< 作業空間プロット点の数 (例：「1S軸～6T軸～7加速度センサ～8力覚センサ～9ツール先端」の9個の要素)
		static constexpr double PLOTXYXZ_XMAX =  1.5;	//!< [m] X軸最大値
		static constexpr double PLOTXYXZ_XMIN = -0.5;	//!< [m] X軸最小値
		static constexpr size_t PLOTXYXZ_XGRID = 4;		//!< X軸グリッドの分割数
		static constexpr char PLOTXYXZ_XLABEL[] = "POSITION X [m]";	//!< X軸ラベル

		//! @brief 作業空間XYプロットの設定
		static constexpr int PLOTXY_LEFT = 1015;		//!< [px] 左位置
		static constexpr int PLOTXY_TOP = 709;			//!< [px] 上位置
		static constexpr int PLOTXY_WIDTH = 355;		//!< [px] 幅
		static constexpr int PLOTXY_HEIGHT = 306;		//!< [px] 高さ
		static constexpr char PLOTXY_YLABEL[] = "POSITION Y [m]";	//!< Y軸ラベル
		static constexpr double PLOTXY_YMAX =  1.0;		//!< [m] Y軸最大値
		static constexpr double PLOTXY_YMIN = -1.0;		//!< [m] Y軸最小値
		static constexpr size_t PLOTXY_YGRID = 4;		//!< Y軸グリッドの分割数
		static constexpr double PLOTXY_VAL_XPOS = -0.4;	//!< 数値表示の左位置
		static constexpr double PLOTXY_VAL_YPOS =  0.9;	//!< 数値表示の上位置
		
		//! @brief 作業空間XZプロットの設定
		static constexpr int PLOTXZ_LEFT = 1370;		//!< [px] 左位置
		static constexpr int PLOTXZ_TOP = 709;			//!< [px] 上位置
		static constexpr int PLOTXZ_WIDTH = 355;		//!< [px] 幅
		static constexpr int PLOTXZ_HEIGHT = 306;		//!< [px] 高さ
		static constexpr char PLOTXZ_ZLABEL[] = "POSITION Z [m]";	//!< Z軸ラベル
		static constexpr double PLOTXZ_ZMAX =  2.0;		//!< [m] Z軸最大値
		static constexpr double PLOTXZ_ZMIN =  0.0;		//!< [m] Z軸最小値
		static constexpr size_t PLOTXZ_ZGRID = 4;		//!< Z軸グリッドの分割数
		static constexpr double PLOTXZ_VAL_XPOS = -0.4;	//!< 数値表示の左位置
		static constexpr double PLOTXZ_VAL_ZPOS =  1.9;	//!< 数値表示の上位置
		
	private:
		ConstParams() = delete;	//!< コンストラクタ使用禁止
		~ConstParams() = delete;//!< デストラクタ使用禁止
		ConstParams(const ConstParams&) = delete;					//!< コピーコンストラクタ使用禁止
		const ConstParams& operator=(const ConstParams&) = delete;	//!< 代入演算子使用禁止
};
}

#endif

