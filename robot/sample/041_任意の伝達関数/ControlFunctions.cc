//! @file ControlFunctions.cc
//! @brief 制御用周期実行関数群クラス
//! @date 2025/01/17
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2025 Yokokura, Yuki
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
#include "ArcsControl.hh"
#include "SquareWave.hh"

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
	constexpr ArcsMat<3,1> N = {5,6,7};		// 伝達関数1の分子の係数ベクトル (5*s^2 + 6*s + 7)
	constexpr ArcsMat<4,1> D = {2,3,4,5};	// 伝達関数1の分母の係数ベクトル (2*s^3 + 3*s^2 + 4*s + 5)
	constexpr ArcsMat<4,1> P = {4,5,6,7};	// 伝達関数2の分子の係数ベクトル (4*s^3 + 5*s^2 + 6*s + 7)
	constexpr ArcsMat<4,1> Q = {2,3,4,5};	// 伝達関数2の分母の係数ベクトル (2*s^3 + 3*s^2 + 4*s + 5)
	constexpr double w = 30;						// [rad/s] 2次遅れ系の周波数
	constexpr double z = 0.2;						// [-] 2次遅れ系の制動係数
	constexpr ArcsMat<1,1> B = {w*w};				// 2次遅れ系の分子係数ベクトル (ω^2)
	constexpr ArcsMat<3,1> A = {1, 2.0*z*w, w*w};	// 2次遅れ系の分母係数ベクトル (s^2 + 2ζωs + ω^2)
	
	// 制御用変数宣言
	static double u = 0, y1 = 0, y2 = 0, y3 = 0;	// 伝達関数の入力信号と出力信号
	static ArcsControl::TransFunc<2,3> G1(N, D, Ts);// 2/3次の伝達関数1の生成
	static ArcsControl::TransFunc<3,3> G2(P, Q, Ts);// 3/3次の伝達関数2の生成(直達項あり)
	static ArcsControl::TransFunc<0,2> G3(B, A, Ts);// 0/2次の伝達関数3の生成(2次標準形)
	
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
		Interface.GetPosition(thm);		// [rad] 位置ベクトルの取得
		Screen.GetOnlineSetVar();		// オンライン設定変数の読み込み
		
		// ここに制御アルゴリズムを記述する
		u = SquareWave(0.5, M_PI, t, 1);// 方形波入力生成
		y1 = G1.GetResponse(u);			// 伝達関数1の出力応答
		y2 = G2.GetResponse(u);			// 伝達関数2の出力応答
		y3 = G3.GetResponse(u);			// 伝達関数3の出力応答

		Interface.SetCurrent(iqref);	// [A] 電流指令ベクトルの出力
		Screen.SetVarIndicator(u, y1, y2, y3, 0, 0, 0, 0, 0, 0);	// 任意変数インジケータ(変数0, ..., 変数9)
		Graph.SetTime(Tact, t);									// [s] グラフ描画用の周期と時刻のセット
		Graph.SetVars(0, u, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット0 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(1, y1, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット1 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(2, y2, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット2 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(3, y3, 0, 0, 0, 0, 0, 0, 0);	// グラフプロット3 (グラフ番号, 変数0, ..., 変数7)
		UsrGraph.SetVars(0, 0);						// ユーザカスタムプロット（例）
		Memory.SetData(Tact, t, u, y1, y2, y3, 0, 0, 0, 0, 0);	// CSVデータ保存変数 (周期, A列, B列, ..., J列)
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

