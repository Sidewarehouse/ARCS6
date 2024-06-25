//! @file ControlFunctions.cc
//! @brief 制御用周期実行関数群クラス
//! @date 2024/06/25
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// 基本のインクルードファイル
#include <cmath>
#include "ARCSprint.hh"
#include "ARCSassert.hh"
#include "ARCSmemory.hh"
#include "ARCSscrparams.hh"
#include "ARCSgraphics.hh"
#include "ControlFunctions.hh"
#include "InterfaceFunctions.hh"

// 追加のARCSライブラリをここに記述
#include "ArcsMatrix.hh"
#include "ZFA01.hh"

using namespace ARCS;

//! @brief スレッド間通信用グローバル変数の無名名前空間
namespace {
	// スレッド間で共有したい変数をここに記述
	ArcsMat<EquipParams::ACTUATOR_NUM, 1> thm;	//!< [rad]  位置ベクトル
	ArcsMat<EquipParams::ACTUATOR_NUM, 1> iqref;//!< [A,Nm] 電流指令,トルク指令ベクトル
}

//! @brief 制御用周期実行関数1
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction1(const double t, const double Tact, const double Tcmp){
	// 制御用定数設定
	[[maybe_unused]] constexpr double Ts = ConstParams::SAMPLING_TIME[0]*1e-9;	// [s]	制御周期
	
	// 制御用変数宣言
	static ArcsMat<ZFA01::N_AX, 1> thl;	// [rad] 負荷側位置ベクトル
	static ArcsMat<ZFA01::N_AX, 1> wm;	// [rad/s] モータ側速度ベクトル
	static ArcsMat<ZFA01::N_AX, 1> wl;	// [rad/s] 負荷側速度ベクトル
	static ArcsMat<ZFA01::N_AX, 1> taus;// [Nm] ねじれトルクベクトル
	
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
		Interface.GetPositionAndVelocity(thm, wm);	// [rad,rad/s] 位置と速度ベクトルの取得
		Interface.GetTorque(taus);					// [Nm] ねじれトルクベクトル応答の取得
		Screen.GetOnlineSetVar();					// オンライン設定変数の読み込み
		
		// ここに制御アルゴリズムを記述する
		thl = thm % ZFA01::Rg;			// [rad] 負荷側位置(簡易的換算)
		wl = wm % ZFA01::Rg;			// [rad/s] 負荷側速度(簡易的換算)
		iqref[1] = 0;					// [A] 1S電流指令
		iqref[2] = 0;					// [A] 2L電流指令
		iqref[3] = 0;					// [A] 3U電流指令
		
		Interface.SetCurrent(iqref);	// [A] 電流指令ベクトルの出力
		Screen.SetVarIndicator(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);	// 任意変数インジケータ(変数0, ..., 変数9)
		Graph.SetTime(Tact, t);									// [s] グラフ描画用の周期と時刻のセット
		Graph.SetVars(0, iqref[1], iqref[2], iqref[3], 0, 0, 0, 0, 0);	// グラフプロット0 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(1, thm[1], thm[2], thm[3], 0, 0, 0, 0, 0);		// グラフプロット1 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(2, wm[1], wm[2], wm[3], 0, 0, 0, 0, 0);			// グラフプロット2 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(3, taus[1], taus[2], taus[3], 0, 0, 0, 0, 0);		// グラフプロット3 (グラフ番号, 変数0, ..., 変数7)
		//Graph.SetWorkspace(0, 0, 0);									// 作業空間グラフプロット
		UsrGraph.SetVars(0, 0);											// ユーザカスタムプロット（例）
		Memory.SetData(Tact, t, iqref[1], taus[1], 0, 0, 0, 0, 0, 0, 0);// CSVデータ保存変数 (周期, A列, B列, ..., J列)
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
bool ControlFunctions::ControlFunction2(const double t, const double Tact, const double Tcmp){
	// 制御用定数宣言
	[[maybe_unused]] constexpr double Ts = ConstParams::SAMPLING_TIME[1]*1e-9;	// [s]	制御周期
	
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
bool ControlFunctions::ControlFunction3(const double t, const double Tact, const double Tcmp){
	// 制御用定数宣言
	[[maybe_unused]] constexpr double Ts = ConstParams::SAMPLING_TIME[2]*1e-9;	// [s]	制御周期
	
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
	Screen.SetNetworkLink(NetworkLink);			// ネットワークリンクフラグを書き込む
	Screen.SetInitializing(Initializing);		// ロボット初期化フラグを書き込む
	Screen.SetCurrentAndPosition(iqref, thm);	// 電流指令と位置を書き込む
}

