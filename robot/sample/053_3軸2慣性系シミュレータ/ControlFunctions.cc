//! @file ControlFunctions.cc
//! @brief 制御用周期実行関数群クラス
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

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
#include "TwoInertiaSimulators.hh"
#include "StepWave.hh"

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
	[[maybe_unused]] constexpr double Ts = ConstParams::SAMPLING_TIME[0]*1e-9;	// [s]	制御周期
	constexpr size_t N = 3;		// [-]		軸数
	constexpr double Kt = 0.49;		// [Nm/A]	トルク定数
	constexpr double Jm = 1.10e-4;	// [kgm^2] モータ側慣性
	constexpr double Dm = 3.76e-4;	// [Nm s/rad]モータ側粘性
	constexpr double Jl = 0.133;	// [kgm^2]	負荷側慣性
	constexpr double Dl = 0.001;	// [Nm s/rad]負荷側粘性
	constexpr double Ks = 1.08e4;  	// [Nm/rad]	2慣性間の剛性
	constexpr double Rg = 50;		// [-]		減速比
	constexpr TwoInertiaParamDef PlantParams = {Kt, Jm, Dm, Jl, Dl, Ks, Rg};	// 2慣性系パラメータ構造体
	constexpr std::array<TwoInertiaParamDef, N> PlantParamsArray = {PlantParams, PlantParams, PlantParams};	// 2慣性系パラメータ構造体の配列
	
	// 制御用変数宣言
	static Matrix<1,N> iq, tl;		// 電流と負荷トルク
	static Matrix<1,N> wl, ths, wm;	// 負荷側速度，ねじれ角，モータ側速度
	static Matrix<1,N> thl, thm;	// 負荷側位置，モータ側位置
	static TwoInertiaSimulators<N> Plant(PlantParamsArray, Ts);	// 2慣性系シミュレータ
	
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
		
		// 3軸分の2慣性系シミュレータの例(ただし全軸独立で干渉は無し)
		iq[1] = StepWave(0.1, t);
		iq[2] = StepWave(0.2, t);
		iq[3] = StepWave(0.3, t);
		Plant.SetCurrentAndLoadTorque(iq, tl);
		Plant.GetResponses(wl, ths, wm);
		Plant.GetLoadPosition(thl);
		Plant.GetMotorPosition(thm);
		
		Interface.SetCurrent(CurrentRef);	// [A] 電流指令の出力
		Screen.SetVarIndicator(wl[1], wl[2], wl[3], thm[1], thm[2], thm[3], 0, 0, 0, 0);	// 任意変数インジケータ(変数0, ..., 変数9)
		Graph.SetTime(Tact, t);									// [s] グラフ描画用の周期と時刻のセット
		Graph.SetVars(0, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット0 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(1, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット1 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(2, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット2 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(3, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット3 (グラフ番号, 変数0, ..., 変数7)
		Memory.SetData(Tact, t, wl[1], wl[2], wl[3], thm[1], thm[2], thm[3], 0, 0, 0);		// CSVデータ保存変数 (周期, A列, B列, ..., J列)
		// リアルタイム制御ここまで
	}
	if(CmdFlag == CTRL_EXIT){
		// 終了処理モード (ここは制御終了時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
		Interface.SetZeroCurrent();	// 電流指令を零に設定
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
	[[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[1]*1e-9;	// [s]	制御周期
	
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
	[[maybe_unused]] const double Ts = ConstParams::SAMPLING_TIME[2]*1e-9;	// [s]	制御周期
	
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

