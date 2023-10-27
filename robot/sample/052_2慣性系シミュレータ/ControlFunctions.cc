//! @file ControlFunctions.cc
//! @brief 制御用周期実行関数群クラス
//! @date 2023/10/26
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// 基本のインクルードファイル
#include <unistd.h>
#include <cmath>
#include <cfloat>
#include <tuple>
#include "ControlFunctions.hh"
#include "ARCSprint.hh"
#include "ARCSassert.hh"
#include "ScreenParams.hh"
#include "InterfaceFunctions.hh"
#include "GraphPlot.hh"
#include "DataMemory.hh"

// 追加のARCSライブラリをここに記述
#include "Matrix.hh"
#include "TwoInertiaParamDef.hh"
#include "TwoInertiaSimulator.hh"
#include "TwoInertiaStateObsrv.hh"
#include "SquareWave.hh"

using namespace ARCS;

//! @brief スレッド間通信用グローバル変数の無名名前空間
namespace {
	// スレッド間で共有したい変数をここに記述
	std::array<double, ConstParams::ACTUATOR_NUM> PositionRes = {0};	//!< [rad] 位置応答
	std::array<double, ConstParams::ACTUATOR_NUM> CurrentRef = {0};		//!< [Nm]  電流指令
}

//! @brief 制御用周期実行関数1
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction1(double t, double Tact, double Tcmp){
	// 制御用定数設定
	[[maybe_unused]] constexpr double Ts = ConstParams::SAMPLING_TIME[0]*1e-9;	// [s]		制御周期
	constexpr double Jl = 0.1;		// [kgm^2]		負荷側慣性
	constexpr double Dl = 1e-12;	// [Nm/(rad/s)]	負荷側粘性
	constexpr double Ds = 1; 	 	// [Nm/(rad/s)]	ねじれ粘性
	constexpr double Ks = 10000;  	// [Nm/rad]		ねじれ剛性
	constexpr double Jm = 1e-5;		// [kgm^2]		モータ側慣性
	constexpr double Dm = 1e-4;		// [Nm/(rad/s)]	モータ側粘性
	constexpr double Rg = 100;		// [-]			減速比
	constexpr double Kt = 0.1;		// [Nm/A]		トルク定数
	constexpr TwoInertiaParams PlantParams = {Jl, Dl, Ds, Ks, Jm, Dm, Rg, Kt};	// 2慣性系パラメータ構造体
	
	// 制御用変数宣言
	static double iq, tdis, wl, ths, wm, taus;			// 電流，外乱, 負荷側速度，ねじれ角，モータ側速度，ねじれトルク
	static TwoInertiaSimulator Plant(PlantParams, Ts);	// 2慣性共振系シミュレータ
	
	if(CmdFlag == CTRL_INIT){
		// 初期化モード (ここは制御開始時/再開時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
		Initializing = true;		// 初期化中ランプ点灯
		Screen.InitOnlineSetVar();	// オンライン設定変数の初期値の設定
		Interface.ServoON();		// サーボON指令の送出
		Initializing = false;		// 初期化中ランプ消灯
	}
	if(CmdFlag == CTRL_LOOP){
		// 周期モード (ここは制御周期 SAMPLING_TIME[0] 毎に呼び出される(リアルタイム空間なので処理は制御周期内に収めること))
		// リアルタイム制御ここから
		Interface.GetPosition(PositionRes);	// [rad] 位置応答の取得
		Screen.GetOnlineSetVar();			// オンライン設定変数の読み込み
		
		PositionRes[0] = Plant.GetMotorPosition();		// [rad] 2慣性系シミュレータからモータ位置を取得
		std::tie(wl, ths, wm) = Plant.GetResponses();	// [rad/s,rad,rad/s] 2慣性系シミュレータから負荷側速度、ねじれ角、モータ側速度を取得
		taus = Plant.GetTorsionTorque();				// [Nm]  2慣性系シミュレータからねじれトルクを取得
		
		iq   = 0.5*SquareWave(0.1, 0, t, 0.5);			// [A]  q軸電流の方形波生成
		tdis =   2*SquareWave(0.2, 0, t, 2.0);			// [Nm] 負荷側外乱トルクの方形波生成
		
		Plant.SetCurrentAndLoadTorque(iq, tdis);		// [A,Nm] 2慣性系シミュレータへの入力と状態更新
		CurrentRef[0] = iq;
		
		Interface.SetTorque(CurrentRef);	// [A] 電流指令の出力
		Screen.SetVarIndicator(wl, ths, wm, taus, 0, 0, 0, 0, 0, 0);	// 任意変数インジケータ(変数0, ..., 変数9)
		Graph.SetTime(Tact, t);				// [s] グラフ描画用の周期と時刻のセット
		Graph.SetVars(0, iq, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット0 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(1, wl, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット1 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(2, ths, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット2 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(3, wm, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット3 (グラフ番号, 変数0, ..., 変数7)
		Memory.SetData(Tact, t, iq, tdis, wl, ths, wm, taus, 0, 0, 0);	// CSVデータ保存変数 (周期, A列, B列, ..., J列)
		// リアルタイム制御ここまで
	}
	if(CmdFlag == CTRL_EXIT){
		// 終了処理モード (ここは制御終了時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
		Interface.SetZeroCurrent();// 電流指令を零に設定
		Interface.ServoOFF();		// サーボOFF信号の送出
	}
	return true;	// クロックオーバーライドフラグ(falseにすると次の周期時刻を待たずにスレッドが即刻動作する)
}

//! @brief 制御用周期実行関数2
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction2(double t, double Tact, double Tcmp){
	// 制御用定数宣言
	[[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[1]*1e-9;	// [s]		制御周期
	
	// 制御用変数宣言
	
	// 制御器等々の宣言
	
	if(CmdFlag == CTRL_INIT){
		// 初期化モード (ここは制御開始時/再開時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	if(CmdFlag == CTRL_LOOP){
		// 周期モード (ここは制御周期 SAMPLING_TIME[1] 毎に呼び出される(リアルタイム空間なので処理は制御周期内に収めること))
		// リアルタイム制御ここから
		
		// リアルタイム制御ここまで
	}
	if(CmdFlag == CTRL_EXIT){
		// 終了処理モード (ここは制御終了時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	return true;	// クロックオーバーライドフラグ(falseにすると次の周期時刻を待たずにスレッドが即刻動作する)
}

//! @brief 制御用周期実行関数3
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction3(double t, double Tact, double Tcmp){
	// 制御用定数宣言
	
	// 制御用変数宣言
	
	if(CmdFlag == CTRL_INIT){
		// 初期化モード (ここは制御開始時/再開時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	if(CmdFlag == CTRL_LOOP){
		// 周期モード (ここは制御周期 SAMPLING_TIME[2] 毎に呼び出される(リアルタイム空間なので処理は制御周期内に収めること))
		// リアルタイム制御ここから
		
		// リアルタイム制御ここまで
	}
	if(CmdFlag == CTRL_EXIT){
		// 終了処理モード (ここは制御終了時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
	}
	return true;	// クロックオーバーライドフラグ(falseにすると次の周期時刻を待たずにスレッドが即刻動作する)
}

//! @brief 制御用変数値を更新する関数
void ControlFunctions::UpdateControlValue(void){
	// ARCS画面パラメータに値を書き込む
	Screen.SetNetworkLink(NetworkLink);						// ネットワークリンクフラグを書き込む
	Screen.SetInitializing(Initializing);					// ロボット初期化フラグを書き込む
	Screen.SetCurrentAndPosition(CurrentRef, PositionRes);	// 電流指令と位置応答を書き込む
}

