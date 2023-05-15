//! @file ConstParams.cc
//! @brief 定数値格納用クラス
//!        ARCSに必要な定数値を格納します。
//! @date 2022/01/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "ConstParams.hh"

// 注意！： const std::string は ConstParams.cc で定義すること。それ以外は constexpr として ConstParams.hh で定義。

using namespace ARCS;

// タイトルに表示させる制御系の名前(識別用に好きな名前を入力)
const std::string ConstParams::CTRLNAME("<TITLE: ZFA01 VELOCITY CONTROL 2L >");		//!< (画面に入る文字数以内)

// ARCS改訂番号(ARCS本体側システムコード改変時にちゃんと変えること)
const std::string ConstParams::ARCS_REVISION("REV.ZFA122011216");	//!< (16文字以内)

// イベントログの設定
const std::string ConstParams::EVENTLOG_NAME("EventLog.txt");		//!< イベントログファイル名

// 実験データCSVファイルの設定
const std::string ConstParams::DATA_NAME("DATA.csv");				//!< CSVファイル名

// SCHED_FIFOリアルタイムスレッドの設定
constexpr std::array<unsigned long, ConstParams::THREAD_MAX> ConstParams::SAMPLING_TIME;	//!< 制御周期の設定
constexpr std::array<unsigned int, ConstParams::THREAD_MAX> ConstParams::CPUCORE_NUMBER;	//!< 使用CPUコアの設定

// 実験機アクチュエータの設定
constexpr std::array<ConstParams::ActType, ConstParams::ACTUATOR_MAX> ConstParams::ACT_TYPE;		//!< 実験機アクチュエータの種類の設定（リニアモータか回転モータかの設定）
constexpr std::array<ConstParams::ActRefUnit, ConstParams::ACTUATOR_MAX> ConstParams::ACT_REFUNIT;	//!< 実験機アクチュエータの指令単位の設定（電流なのか推力なのかトルクなのかの設定）
constexpr std::array<double, ConstParams::ACTUATOR_MAX> ConstParams::ACT_FORCE_TORQUE_CONST;//!< トルク/推力定数の設定
constexpr std::array<double, ConstParams::ACTUATOR_MAX> ConstParams::ACT_RATED_CURRENT;		//!< 定格電流値の設定
constexpr std::array<double, ConstParams::ACTUATOR_MAX> ConstParams::ACT_MAX_CURRENT;		//!< 瞬時最大許容電流値の設定
constexpr std::array<double, ConstParams::ACTUATOR_MAX> ConstParams::ACT_RATED_TORQUE;		//!< 定格トルクの設定
constexpr std::array<double, ConstParams::ACTUATOR_MAX> ConstParams::ACT_MAX_TORQUE;		//!< 瞬時最大トルクの設定
constexpr std::array<double, ConstParams::ACTUATOR_MAX> ConstParams::ACT_INITPOS;			//!< 初期位置の設定

// 任意変数値表示の設定
const std::string ConstParams::INDICVARS_FORMS[INDICVARS_MAX] = {
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
	"% 13.4f",	// 任意に表示したい変数値の表示形式(printfの書式と同一)
};

// 時系列グラフプロットの共通設定
const std::string ConstParams::PLOT_FRAMEBUFF("/dev/fb0");			//!< フレームバッファ ファイルデスクリプタ
const std::string ConstParams::PLOT_PNGFILENAME("Screenshot.png");	//!< スクリーンショットのPNGファイル名
const std::string ConstParams::PLOT_TFORMAT("%3.1f");				//!< 横軸書式
const std::string ConstParams::PLOT_TLABEL("Time [s]");				//!< 横軸ラベル
constexpr std::array<bool, ConstParams::PLOT_MAX> ConstParams::PLOT_VISIBLE;			//!< グラフ描画の有効/無効設定
constexpr std::array<FGcolors, ConstParams::PLOT_VAR_MAX> ConstParams::PLOT_VAR_COLORS;	//!< 線の色
constexpr std::array<unsigned int, ConstParams::PLOT_MAX> ConstParams::PLOT_VAR_NUM;	//!< プロットする変数の数 (≦PLOT_VAR_MAX)
constexpr std::array<double, ConstParams::PLOT_MAX> ConstParams::PLOT_FMAX;				//!< 縦軸最大値
constexpr std::array<double, ConstParams::PLOT_MAX> ConstParams::PLOT_FMIN;				//!< 縦軸最小値
constexpr std::array<unsigned int, ConstParams::PLOT_MAX> ConstParams::PLOT_FGRID_NUM;	//!< 縦軸グリッドの分割数
constexpr std::array<int, ConstParams::PLOT_MAX> ConstParams::PLOT_LEFT;				//!< [px] 左位置
constexpr std::array<int, ConstParams::PLOT_MAX> ConstParams::PLOT_TOP;					//!< [px] 上位置
constexpr std::array<int, ConstParams::PLOT_MAX> ConstParams::PLOT_WIDTH;				//!< [px] 幅
constexpr std::array<int, ConstParams::PLOT_MAX> ConstParams::PLOT_HEIGHT;				//!< [px] 高さ
constexpr std::array<
	std::array<CuiPlotTypes, ConstParams::PLOT_VAR_MAX>, ConstParams::PLOT_MAX
> ConstParams::PLOT_TYPE;																//!< プロットの種類

// 時系列グラフプロットの各部設定
//!< @brief 縦軸ラベルの設定
const std::array<std::string, ConstParams::PLOT_MAX> ConstParams::PLOT_FLABEL = {
	"CURRENT REF [A]",		// プロット0
	"POSITION [rad]",		// プロット1
	"VELOCITY [rad/s]",		// プロット2
	"TORSION TORQUE [Nm]",	// プロット3
	"WORKSPACE VELOCITY [m/s]",	// プロット4
	"WORKSPACE FORCE [N]",		// プロット5
	"---------- [-]",	// プロット6
	"---------- [-]",	// プロット7
	"---------- [-]",	// プロット8
	"---------- [-]",	// プロット9
	"---------- [-]",	// プロット10
	"---------- [-]",	// プロット11
	"---------- [-]",	// プロット12
	"---------- [-]",	// プロット13
	"---------- [-]",	// プロット14
	"---------- [-]",	// プロット15
};

//! @brief 時系列グラフプロットの縦軸書式
const std::array<std::string, ConstParams::PLOT_MAX> ConstParams::PLOT_FFORMAT = {
	"%6.1f",	// プロット0
	"%6.0f",	// プロット1
	"%6.0f",	// プロット2
	"%6.0f",	// プロット3
	"%6.1f",	// プロット4
	"%6.0f",	// プロット5
	"%6.1f",	// プロット6
	"%6.1f",	// プロット7
	"%6.1f",	// プロット8
	"%6.1f",	// プロット9
	"%6.1f",	// プロット10
	"%6.1f",	// プロット11
	"%6.1f",	// プロット12
	"%6.1f",	// プロット13
	"%6.1f",	// プロット14
	"%6.1f",	// プロット15
};

//! @brief 時系列グラフプロットの変数の名前
const std::array<
	std::array<std::string, ConstParams::PLOT_VAR_MAX>,
	ConstParams::PLOT_MAX
> ConstParams::PLOT_VAR_NAMES = {{
	{"1S", "2L", "3U", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット0
	{"1S", "2L", "3U", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット1
	{"1S", "2L", "3U", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット2
	{"1S", "2L", "3U", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット3
	{"vx", "vy", "vz", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット4
	{"fx", "fy", "fz", "VAR-03", "VAR-04", "VAR-05", "VAR-06", "VAR-07",},	// プロット5
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

// 作業空間XYプロットの設定
const std::string ConstParams::PLOTXY_XLABEL("POSITION X [m]");		//!< X軸ラベル
const std::string ConstParams::PLOTXY_YLABEL("POSITION Y [m]");		//!< Y軸ラベル

// 作業空間XZプロットの設定
const std::string ConstParams::PLOTXZ_XLABEL("POSITION X [m]");		//!< X軸ラベル
const std::string ConstParams::PLOTXZ_ZLABEL("POSITION Z [m]");		//!< Z軸ラベル

// ユーザプロットの設定
const std::string ConstParams::PLOTUS_XLABEL("X AXIS [-]");	//!< X軸ラベル
const std::string ConstParams::PLOTUS_YLABEL("Y AXIS [-]");	//!< Y軸ラベル
