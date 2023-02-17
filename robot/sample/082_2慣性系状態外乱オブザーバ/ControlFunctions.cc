//! @file ControlFunctions.cc
//! @brief 制御用周期実行関数群クラス
//! @date 2020/04/06
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
#include "TwoInertiaParamDef.hh"
#include "TwoInertiaSimulator.hh"
#include "TwoInertiaStateDistObsrv.hh"
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
	constexpr double Kt = 0.49;		// [Nm/A]	トルク定数
	constexpr double Jm = 1.10e-4;	// [kgm^2] モータ側慣性
	constexpr double Dm = 3.76e-4;	// [Nm s/rad]モータ側粘性
	constexpr double Jl = 0.133;	// [kgm^2]	負荷側慣性
	constexpr double Dl = 0.1;		// [Nm s/rad]負荷側粘性
	constexpr double Ks = 1.08e4;  	// [Nm/rad]	2慣性間の剛性
	constexpr double Rg = 50;		// [-]		減速比
	constexpr TwoInertiaParamDef PlantParams = {Kt, Jm, Dm, Jl, Dl, Ks, Rg};	// 2慣性系パラメータ構造体
	
	// 制御用変数宣言
	static double iq, wl, ths, wm, tl;	// 電流，速度，位置，負荷側外乱
	static Matrix<1,4> xhat;			// 推定状態ベクトル
	static Matrix<1,3> xhat2;			// 推定状態ベクトル
	static double wlh, thsh, wmh, tlh;	// 推定値
	static TwoInertiaSimulator Plant(PlantParams, Ts);			// 2慣性共振系シミュレータ
	static TwoInertiaStateDistObsrv SDOB(PlantParams, 300, Ts);	// 状態外乱オブザーバ
	
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
		
		iq = SquareWave(1, 0, t);							// [A] q軸電流の方形波
		tl = SquareWave(2, M_PI, t);						// [Nm]負荷側のステップ外乱
		std::tie(wl, ths, wm) = Plant.GetResponses(iq, tl);	// [rad/s,rad] 2慣性系シミュレータ
		//xhat = SDOB.GetEstimatedVect(iq, wm);							// [rad/s,rad] 状態外乱オブザーバ(状態ベクトルが欲しい場合)
		//std::tie(wlh, thsh, wmh, tlh) = SDOB.GetEstimatedVars(iq, wm);// [rad/s,rad] 状態外乱オブザーバ(状態変数それぞれが欲しい場合)
		SDOB.GetEstimatedVect(iq, wm, xhat2, tlh);			// [rad/s,rad] 状態外乱オブザーバ(状態ベクトルと負荷側外乱トルクを分けて欲しい場合)
		wlh  = xhat2[1];	// [rad/s] ωl推定値
		thsh = xhat2[2];	// [rad] θs推定値
		wmh  = xhat2[3];	// [rad/s] ωm推定値
		
		Interface.SetTorque(CurrentRef);	// [Nm] トルク指令の出力
		Screen.SetVarIndicator(wl, wlh, ths, thsh, wm, wmh, tl, tlh, 0, 0);	// 任意変数インジケータ(変数0, ..., 変数9)
		Graph.SetTime(Tact, t);				// [s] グラフ描画用の周期と時刻のセット
		Graph.SetVars(0, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット0 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(1, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット1 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(2, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット2 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(3, 0, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット3 (グラフ番号, 変数0, ..., 変数7)
		Memory.SetData(Tact, t, wl, wlh, ths, thsh, wm, wmh, tl, tlh, 0);		// CSVデータ保存変数 (周期, A列, B列, ..., J列)
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

