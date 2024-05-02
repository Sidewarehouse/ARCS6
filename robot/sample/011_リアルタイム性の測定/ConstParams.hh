//! @file ConstParams.hh
//! @brief 定数値格納用クラス
//!        ARCSに必要な定数値を格納します。
//! @date 2024/05/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef CONSTPARAMS
#define CONSTPARAMS

#include <cmath>
#include "SFthread.hh"
#include "FrameGraphics.hh"
#include "CuiPlot.hh"

// 注意！： const std::string は ConstParams.cc で定義すること。それ以外は constexpr として ConstParams.hh で定義。

namespace ARCS {	// ARCS名前空間
	//! @brief 定数値格納用クラス
	class ConstParams {
		public:
			// タイトルに表示させる制御系の名前(識別用に好きな名前を入力)
			static const std::string CTRLNAME;		//!< (60文字以内)
			
			// 画面サイズの設定 (モニタ解像度に合うように設定すること)
			// 1024×600(WSVGA) の場合に下記をアンコメントすること
			//static constexpr int SCR_VERTICAL_MAX = 36;		//!< [文字] 画面の最大高さ文字数
			//static constexpr int SCR_HORIZONTAL_MAX = 127;	//!< [文字] 画面の最大幅文字数
			// 1024×768(XGA) の場合に下記をアンコメントすること
			static constexpr int SCR_VERTICAL_MAX = 47;			//!< [文字] 画面の最大高さ文字数
			static constexpr int SCR_HORIZONTAL_MAX = 127;		//!< [文字] 画面の最大幅文字数
			// 1280×1024(SXGA) の場合に下記をアンコメントすること
			//static constexpr int SCR_VERTICAL_MAX = 63;		//!< [文字] 画面の最大高さ文字数
			//static constexpr int SCR_HORIZONTAL_MAX = 159;	//!< [文字] 画面の最大幅文字数
			// 1920×1080(Full HD) の場合に下記をアンコメントすること
			//static constexpr int SCR_VERTICAL_MAX = 66;		//!< [文字] 画面の最大高さ文字数
			//static constexpr int SCR_HORIZONTAL_MAX = 239;	//!< [文字] 画面の最大幅文字数
			// それ以外の解像度の場合は各自で値を探すこと
			
			// 実験データCSVファイルの設定
			static const std::string DATA_NAME;				//!< CSVファイル名
			static constexpr double DATA_START =  0;		//!< [s] 保存開始時刻
			static constexpr double DATA_END   = 10;		//!< [s] 保存終了時刻
			static constexpr double DATA_RESO  = 0.001;		//!< [s] データの時間分解能
			static constexpr unsigned int DATA_NUM  =  10;	//!< [-] 保存する変数の数
			
			// SCHED_FIFOリアルタイムスレッドの設定
			static constexpr unsigned int THREAD_MAX = 3;	//!< スレッド最大数（これ変えても ControlFunctions.cc は追随しないので注意）
			static constexpr unsigned int THREAD_NUM = 1;	//!< 動作させるスレッドの数 (最大数は THREAD_NUM_MAX 個まで)
			static constexpr SFalgorithm THREAD_TYPE = SFalgorithm::INSERT_ZEROSLEEP;	//!< リアルタイムアルゴリズムの選択
			// 上記を INSERT_ZEROSLEEP にすると安定性が増すがリアルタイム性は落ちる。遅い処理系の場合に推奨。
			// WITHOUT_ZEROSLEEP にするとリアルタイム性が向上するが，一時的に操作不能になる可能性が残る。高速な処理系の場合に選択可。
			// 下記はカーネルパラメータの設定
			// NO_SETTINGS もしくは CFS_DISABLED と PREEMPT_DYNFULL が併用可。詳細はSFthreadクラスのコメント欄を参照のこと。
			//static constexpr SFkernelparam THREAD_KP = SFkernelparam::CFS_DISABLED;		//!< CFSをリアルタイム用に設定
			static constexpr SFkernelparam THREAD_KP = static_cast<SFkernelparam>(
				static_cast<uint8_t>(SFkernelparam::CFS_DISABLED) | static_cast<uint8_t>(SFkernelparam::PREEMPT_DYNFULL)
			);	//!< CFSとPREEMPTの設定を併用する場合の例
			
			//! @brief 制御周期の設定
			static constexpr std::array<unsigned long, THREAD_MAX> SAMPLING_TIME = {
			//   s  m  u  n	制御周期は Ts[0] ≦ Ts[1] ≦ … ≦ Ts[THREAD_MAX] になるようにすること
				     100000,	// [ns] 制御用周期実行関数1 (スレッド1) 制御周期
				    1000000,	// [ns] 制御用周期実行関数2 (スレッド2) 制御周期
				    1000000,	// [ns] 制御用周期実行関数3 (スレッド3) 制御周期
			};
			
			//! @brief 使用CPUコアの設定
			static constexpr std::array<unsigned int, THREAD_MAX> CPUCORE_NUMBER = {
				    3,	// [-] 制御用周期実行関数1 (スレッド1) 使用するCPUコア番号
				    2,	// [-] 制御用周期実行関数2 (スレッド2) 使用するCPUコア番号
				    1,	// [-] 制御用周期実行関数3 (スレッド3) 使用するCPUコア番号
			};
			
			// 実験機アクチュエータの設定
			static constexpr unsigned int ACTUATOR_MAX = 16;	//!< ARCSが対応しているアクチュエータの最大数
			static constexpr unsigned int ACTUATOR_NUM = 1;		//!< 実験装置のアクチュエータの総数
			
			//! @brief アクチュエータタイプの定義
			enum ActType {
				LINEAR_MOTOR,	//!< リニアモータ
				ROTARY_MOTOR	//!< 回転モータ
			};
			
			//! @brief 実験機アクチュエータの種類の設定（リニアモータか回転モータかの設定）
			static constexpr std::array<ActType, ACTUATOR_MAX> ACT_TYPE = {
				ROTARY_MOTOR,	//  1番 アクチュエータ
				ROTARY_MOTOR,	//  2番 アクチュエータ
				ROTARY_MOTOR,	//  3番 アクチュエータ
				ROTARY_MOTOR,	//  4番 アクチュエータ
				ROTARY_MOTOR,	//  5番 アクチュエータ
				ROTARY_MOTOR,	//  6番 アクチュエータ
				ROTARY_MOTOR,	//  7番 アクチュエータ
				ROTARY_MOTOR,	//  8番 アクチュエータ
				ROTARY_MOTOR,	//  9番 アクチュエータ
				ROTARY_MOTOR,	// 10番 アクチュエータ
				ROTARY_MOTOR,	// 11番 アクチュエータ
				ROTARY_MOTOR,	// 12番 アクチュエータ
				ROTARY_MOTOR,	// 13番 アクチュエータ
				ROTARY_MOTOR,	// 14番 アクチュエータ
				ROTARY_MOTOR,	// 15番 アクチュエータ
				ROTARY_MOTOR,	// 16番 アクチュエータ
			};
			
			//! @brief アクチュエータ指令単位の定義
			enum ActRefUnit {
				AMPERE,			//!< アンペア単位
				NEWTON,			//!< ニュートン単位
				NEWTON_METER	//!< ニュートンメートル単位
			};
			
			//! @brief 実験機アクチュエータの指令単位の設定（電流なのか推力なのかトルクなのかの設定）
			static constexpr std::array<ActRefUnit, ACTUATOR_MAX> ACT_REFUNIT = {
				AMPERE,	//  1番 アクチュエータ
				AMPERE,	//  2番 アクチュエータ
				AMPERE,	//  3番 アクチュエータ
				AMPERE,	//  4番 アクチュエータ
				AMPERE,	//  5番 アクチュエータ
				AMPERE,	//  6番 アクチュエータ
				AMPERE,	//  7番 アクチュエータ
				AMPERE,	//  8番 アクチュエータ
				AMPERE,	//  9番 アクチュエータ
				AMPERE,	// 10番 アクチュエータ
				AMPERE,	// 11番 アクチュエータ
				AMPERE,	// 12番 アクチュエータ
				AMPERE,	// 13番 アクチュエータ
				AMPERE,	// 14番 アクチュエータ
				AMPERE,	// 15番 アクチュエータ
				AMPERE,	// 16番 アクチュエータ
			};
			
			//! @brief トルク/推力定数の設定
			static constexpr std::array<double, ACTUATOR_MAX> ACT_FORCE_TORQUE_CONST = {
				1,	// [N/A]/[Nm/A]  1番 アクチュエータ
				1,	// [N/A]/[Nm/A]  2番 アクチュエータ
				1,	// [N/A]/[Nm/A]  3番 アクチュエータ
				1,	// [N/A]/[Nm/A]  4番 アクチュエータ
				1,	// [N/A]/[Nm/A]  5番 アクチュエータ
				1,	// [N/A]/[Nm/A]  6番 アクチュエータ
				1,	// [N/A]/[Nm/A]  7番 アクチュエータ
				1,	// [N/A]/[Nm/A]  8番 アクチュエータ
				1,	// [N/A]/[Nm/A]  9番 アクチュエータ
				1,	// [N/A]/[Nm/A] 10番 アクチュエータ
				1,	// [N/A]/[Nm/A] 11番 アクチュエータ
				1,	// [N/A]/[Nm/A] 12番 アクチュエータ
				1,	// [N/A]/[Nm/A] 13番 アクチュエータ
				1,	// [N/A]/[Nm/A] 14番 アクチュエータ
				1,	// [N/A]/[Nm/A] 15番 アクチュエータ
				1,	// [N/A]/[Nm/A] 16番 アクチュエータ
			};
			
			//! @brief 定格電流値の設定
			static constexpr std::array<double, ACTUATOR_MAX> ACT_RATED_CURRENT = {
				1,	// [A]  1番 アクチュエータ
				1,	// [A]  2番 アクチュエータ
				1,	// [A]  3番 アクチュエータ
				1,	// [A]  4番 アクチュエータ
				1,	// [A]  5番 アクチュエータ
				1,	// [A]  6番 アクチュエータ
				1,	// [A]  7番 アクチュエータ
				1,	// [A]  8番 アクチュエータ
				1,	// [A]  9番 アクチュエータ
				1,	// [A] 10番 アクチュエータ
				1,	// [A] 11番 アクチュエータ
				1,	// [A] 12番 アクチュエータ
				1,	// [A] 13番 アクチュエータ
				1,	// [A] 14番 アクチュエータ
				1,	// [A] 15番 アクチュエータ
				1,	// [A] 16番 アクチュエータ
			};
			
			//! @brief 瞬時最大許容電流値の設定
			static constexpr std::array<double, ACTUATOR_MAX> ACT_MAX_CURRENT = {
				3,	// [A]  1番 アクチュエータ
				3,	// [A]  2番 アクチュエータ
				3,	// [A]  3番 アクチュエータ
				3,	// [A]  4番 アクチュエータ
				3,	// [A]  5番 アクチュエータ
				3,	// [A]  6番 アクチュエータ
				3,	// [A]  7番 アクチュエータ
				3,	// [A]  8番 アクチュエータ
				3,	// [A]  9番 アクチュエータ
				3,	// [A] 10番 アクチュエータ
				3,	// [A] 11番 アクチュエータ
				3,	// [A] 12番 アクチュエータ
				3,	// [A] 13番 アクチュエータ
				3,	// [A] 14番 アクチュエータ
				3,	// [A] 15番 アクチュエータ
				3,	// [A] 16番 アクチュエータ
			};
			
			//! @brief 定格トルクの設定
			static constexpr std::array<double, ACTUATOR_MAX> ACT_RATED_TORQUE = {
				1,	// [Nm]  1番 アクチュエータ
				1,	// [Nm]  2番 アクチュエータ
				1,	// [Nm]  3番 アクチュエータ
				1,	// [Nm]  4番 アクチュエータ
				1,	// [Nm]  5番 アクチュエータ
				1,	// [Nm]  6番 アクチュエータ
				1,	// [Nm]  7番 アクチュエータ
				1,	// [Nm]  8番 アクチュエータ
				1,	// [Nm]  9番 アクチュエータ
				1,	// [Nm] 10番 アクチュエータ
				1,	// [Nm] 11番 アクチュエータ
				1,	// [Nm] 12番 アクチュエータ
				1,	// [Nm] 13番 アクチュエータ
				1,	// [Nm] 14番 アクチュエータ
				1,	// [Nm] 15番 アクチュエータ
				1,	// [Nm] 16番 アクチュエータ
			};
			
			//! @brief 瞬時最大トルクの設定
			static constexpr std::array<double, ACTUATOR_MAX> ACT_MAX_TORQUE = {
				3,	// [Nm]  1番 アクチュエータ
				3,	// [Nm]  2番 アクチュエータ
				3,	// [Nm]  3番 アクチュエータ
				3,	// [Nm]  4番 アクチュエータ
				3,	// [Nm]  5番 アクチュエータ
				3,	// [Nm]  6番 アクチュエータ
				3,	// [Nm]  7番 アクチュエータ
				3,	// [Nm]  8番 アクチュエータ
				3,	// [Nm]  9番 アクチュエータ
				3,	// [Nm] 10番 アクチュエータ
				3,	// [Nm] 11番 アクチュエータ
				3,	// [Nm] 12番 アクチュエータ
				3,	// [Nm] 13番 アクチュエータ
				3,	// [Nm] 14番 アクチュエータ
				3,	// [Nm] 15番 アクチュエータ
				3,	// [Nm] 16番 アクチュエータ
			};
			
			//! @brief 初期位置の設定
			static constexpr std::array<double, ACTUATOR_MAX> ACT_INITPOS = {
				0,	// [rad]  1軸 アクチュエータ
				0,	// [rad]  2軸 アクチュエータ
				0,	// [rad]  3軸 アクチュエータ
				0,	// [rad]  4軸 アクチュエータ
			    0,	// [rad]  5軸 アクチュエータ
				0,	// [rad]  6軸 アクチュエータ 
				0,	// [rad]  7番 アクチュエータ
				0,	// [rad]  8番 アクチュエータ
				0,	// [rad]  9番 アクチュエータ
				0,	// [rad] 10番 アクチュエータ
				0,	// [rad] 11番 アクチュエータ
				0,	// [rad] 12番 アクチュエータ
				0,	// [rad] 13番 アクチュエータ
				0,	// [rad] 14番 アクチュエータ
				0,	// [rad] 15番 アクチュエータ
				0,	// [rad] 16番 アクチュエータ
			};
			
			// デバッグプリントとデバッグインジケータの設定
			static constexpr bool DEBUG_PRINT_VISIBLE = false;	//!< デバッグプリント表示の有効/無効設定
			static constexpr bool DEBUG_INDIC_VISIBLE = false;	//!< デバッグインジケータ表示の有効/無効設定
			
			// 任意変数値表示の設定
			static constexpr unsigned int INDICVARS_MAX = 16;	//!< 表示変数最大数 (変更不可)
			static constexpr unsigned int INDICVARS_NUM = 10;	//!< 表示したい変数の数 (最大数 INDICVARS_MAX まで)
			static const std::string INDICVARS_FORMS[];			//!< 任意に表示したい変数値の表示形式 (printfの書式と同一)
			
			// オンライン設定変数の設定
			static constexpr unsigned int ONLINEVARS_MAX = 16;	//!< オンライン設定変数最大数 (変更不可)
			static constexpr unsigned int ONLINEVARS_NUM = 10;	//!< オンライン設定変数の数 (最大数 ONLINEVARS_MAX まで)
			
			// 時系列グラフプロットの共通設定
			static const std::string PLOT_FRAMEBUFF;					//!< フレームバッファ ファイルデスクリプタ
			static const std::string PLOT_PNGFILENAME;					//!< スクリーンショットのPNGファイル名
			static constexpr size_t PLOT_MAX = 16;						//!< [-] グラフプロットの最大数 (変更不可)
			static constexpr size_t PLOT_NUM =  4;						//!< [-] グラフプロットの数
			static constexpr double PLOT_TIMESPAN = 10;					//!< [s] プロットの時間幅
			static constexpr double PLOT_TIMERESO = 0.01;				//!< [s] プロットの時間分解能
			static constexpr size_t PLOT_RINGBUFF = 1024;				//!< [-] プロット用リングバッファの要素数
			static constexpr unsigned int PLOT_TGRID_NUM = 10;			//!< [-] 時間軸グリッドの分割数
			static constexpr unsigned int PLOT_VAR_MAX = 8;				//!< [-] プロット可能な変数の最大数 (変更不可)
			static const std::string PLOT_TFORMAT;						//!< 時間軸書式
			static const std::string PLOT_TLABEL;						//!< 時間軸ラベル
			static const std::array<std::string, PLOT_MAX> PLOT_FLABEL;	//!< 縦軸ラベル
			static const std::array<std::string, PLOT_MAX> PLOT_FFORMAT;//!< 縦軸書式
			static const std::array<
				std::array<std::string, PLOT_VAR_MAX>, PLOT_MAX
			> PLOT_VAR_NAMES;											//!< 変数の名前
			static constexpr FGcolors PLOT_AXIS_COLOR = FGcolors::WHITE;	//!< 軸の色
			static constexpr FGcolors PLOT_GRID_COLOR = FGcolors::GRAY25;	//!< グリッドの色
			static constexpr FGcolors PLOT_BACK_COLOR = FGcolors::BLACK;	//!< 背景色
			static constexpr FGcolors PLOT_TEXT_COLOR = FGcolors::WHITE;	//!< 文字色
			static constexpr FGcolors PLOT_CURS_COLOR = FGcolors::GRAY50;	//!< 時刻カーソルの色
			
			//! @brief 時系列グラフ描画の有効/無効設定
			static constexpr std::array<bool, ConstParams::PLOT_MAX> PLOT_VISIBLE = {
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
				false,	// プロット12
				false,	// プロット13
				false,	// プロット14
				false,	// プロット15
			};
			
			//! @brief 時系列プロットの変数ごとの線の色
			static constexpr std::array<FGcolors, PLOT_VAR_MAX> PLOT_VAR_COLORS = {
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
			static constexpr std::array<unsigned int, PLOT_MAX> PLOT_VAR_NUM = {
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
			static constexpr std::array<double, PLOT_MAX> PLOT_FMAX	= {
				100,	// プロット0
				200,	// プロット1
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
			static constexpr std::array<double, PLOT_MAX> PLOT_FMIN = {
				 0.0,	// プロット0
				 0.0,	// プロット1
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
			static constexpr std::array<unsigned int, PLOT_MAX> PLOT_FGRID_NUM = {
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
			static constexpr std::array<int, ConstParams::PLOT_MAX> PLOT_LEFT = {
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
			static constexpr std::array<int, ConstParams::PLOT_MAX> PLOT_TOP = {
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
			static constexpr std::array<int, ConstParams::PLOT_MAX> PLOT_WIDTH = {
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
			static constexpr std::array<int, ConstParams::PLOT_MAX> PLOT_HEIGHT = {
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
			
			//! @brief 時系列プロットの種類
			//! 下記のプロット方法が使用可能
			//!	PLOT_LINE		線プロット
			//!	PLOT_BOLDLINE 	太線プロット
			//!	PLOT_DOT		点プロット
			//!	PLOT_BOLDDOT	太点プロット
			//!	PLOT_CROSS		十字プロット
			//!	PLOT_STAIRS		階段プロット
			//!	PLOT_BOLDSTAIRS	太線階段プロット
			//!	PLOT_LINEANDDOT	線と点の複合プロット
			static constexpr std::array<std::array<CuiPlotTypes, PLOT_VAR_MAX>, PLOT_MAX> PLOT_TYPE = {{
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
			
			//! @brief 作業空間XYプロットの設定
			static constexpr bool PLOTXY_VISIBLE = false;	//!< プロット可視/不可視設定
			static constexpr int PLOTXY_LEFT = 305;			//!< [px] 左位置
			static constexpr int PLOTXY_TOP = 709;			//!< [px] 上位置
			static constexpr int PLOTXY_WIDTH = 300;		//!< [px] 幅
			static constexpr int PLOTXY_HEIGHT = 270;		//!< [px] 高さ
			static const std::string PLOTXY_XLABEL;			//!< X軸ラベル
			static const std::string PLOTXY_YLABEL;			//!< Y軸ラベル
			static constexpr double PLOTXY_XMAX =  1.5;		//!< [m] X軸最大値
			static constexpr double PLOTXY_XMIN = -0.5;		//!< [m] X軸最小値
			static constexpr double PLOTXY_YMAX =  1.0;		//!< [m] Y軸最大値
			static constexpr double PLOTXY_YMIN = -1.0;		//!< [m] Y軸最小値
			static constexpr unsigned int PLOTXY_XGRID = 4;	//!< X軸グリッドの分割数
			static constexpr unsigned int PLOTXY_YGRID = 4;	//!< Y軸グリッドの分割数
			static constexpr double PLOTXY_VAL_XPOS = -0.4;	//!< 数値表示の左位置
			static constexpr double PLOTXY_VAL_YPOS =  0.9;	//!< 数値表示の上位置
			
			//! @brief 作業空間XZプロットの設定
			static constexpr bool PLOTXZ_VISIBLE = false;	//!< プロット可視/不可視設定
			static constexpr int PLOTXZ_LEFT = 605;			//!< [px] 左位置
			static constexpr int PLOTXZ_TOP = 709;			//!< [px] 上位置
			static constexpr int PLOTXZ_WIDTH = 300;		//!< [px] 幅
			static constexpr int PLOTXZ_HEIGHT = 270;		//!< [px] 高さ
			static const std::string PLOTXZ_XLABEL;			//!< X軸ラベル
			static const std::string PLOTXZ_ZLABEL;			//!< Z軸ラベル
			static constexpr double PLOTXZ_XMAX =  1.5;		//!< [m] X軸最大値
			static constexpr double PLOTXZ_XMIN = -0.5;		//!< [m] X軸最小値
			static constexpr double PLOTXZ_ZMAX =  2.0;		//!< [m] Z軸最大値
			static constexpr double PLOTXZ_ZMIN =  0.0;		//!< [m] Z軸最小値
			static constexpr unsigned int PLOTXZ_XGRID = 4;	//!< X軸グリッドの分割数
			static constexpr unsigned int PLOTXZ_ZGRID = 4;	//!< Z軸グリッドの分割数
			static constexpr double PLOTXZ_VAL_XPOS = -0.4;	//!< 数値表示の左位置
			static constexpr double PLOTXZ_VAL_ZPOS =  1.9;	//!< 数値表示の上位置
			
			//! @brief ユーザプロットの設定
			static constexpr bool PLOTUS_VISIBLE = false;	//!< プロット可視/不可視設定
			static constexpr int PLOTUS_LEFT = 905;			//!< [px] 左位置
			static constexpr int PLOTUS_TOP = 709;			//!< [px] 上位置
			static constexpr int PLOTUS_WIDTH = 300;		//!< [px] 幅
			static constexpr int PLOTUS_HEIGHT = 270;		//!< [px] 高さ
			static const std::string PLOTUS_XLABEL;			//!< X軸ラベル
			static const std::string PLOTUS_YLABEL;			//!< Y軸ラベル
			static constexpr double PLOTUS_XMAX =  10;		//!< [mm] X軸最大値
			static constexpr double PLOTUS_XMIN = -10;		//!< [mm] X軸最小値
			static constexpr double PLOTUS_YMAX =  20;		//!< [mm] Y軸最大値
			static constexpr double PLOTUS_YMIN =   0;		//!< [mm] Y軸最小値
			static constexpr unsigned int PLOTUS_XGRID = 4;	//!< X軸グリッドの分割数
			static constexpr unsigned int PLOTUS_YGRID = 4;	//!< Y軸グリッドの分割数
			
		private:
			ConstParams() = delete;	//!< コンストラクタ使用禁止
			~ConstParams() = delete;//!< デストラクタ使用禁止
			ConstParams(const ConstParams&) = delete;					//!< コピーコンストラクタ使用禁止
			const ConstParams& operator=(const ConstParams&) = delete;	//!< 代入演算子使用禁止
	};
}

#endif

