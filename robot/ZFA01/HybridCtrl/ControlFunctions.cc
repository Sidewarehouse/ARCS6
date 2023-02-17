//! @file ControlFunctions.cc
//! @brief 制御用周期実行関数群クラス
//! @date 2022/03/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

// 基本のインクルードファイル
#include <cmath>
#include "ControlFunctions.hh"
#include "ARCSprint.hh"
#include "ARCSassert.hh"
#include "ScreenParams.hh"
#include "InterfaceFunctions.hh"
#include "GraphPlot.hh"
#include "DataMemory.hh"

// 追加のARCSライブラリをここに記述
#include "Matrix.hh"
#include "../robot/ZFA01/addon/ZFA01.hh"
#include "DisturbanceObsrv.hh"
#include "SshapeGenerators.hh"
#include "PulseWave.hh"
#include "HighPassFilter.hh"

using namespace ARCS;

//! @brief スレッド間通信用グローバル変数の無名名前空間
namespace {
	// スレッド間で共有したい変数をここに記述
	std::array<double, ConstParams::ACTUATOR_NUM> PositionRes = {0};//!< [rad] 位置応答
	std::array<double, ConstParams::ACTUATOR_NUM> VelocityRes = {0};//!< [rad/s] 速度応答
	std::array<double, ConstParams::ACTUATOR_NUM> TorqueRes = {0};	//!< [Nm] トルク応答
	std::array<double, ConstParams::ACTUATOR_NUM> CurrentRef = {0};	//!< [Nm]  電流指令
	static std::array<double, 3> ForceRes = {0};					//!< [N,Nm] 6軸力覚応答の並進分
}

//! @brief 制御用周期実行関数1
//! @param[in]	t		時刻 [s]
//! @param[in]	Tact	計測周期 [s]
//! @param[in]	Tcmp	消費時間 [s]
//! @return		クロックオーバーライドフラグ (true = リアルタイムループ, false = 非リアルタイムループ)
bool ControlFunctions::ControlFunction1(double t, double Tact, double Tcmp){
	// 制御用定数設定
	[[maybe_unused]] constexpr double Ts = ConstParams::SAMPLING_TIME[0]*1e-9;	// [s]	制御周期
	static constexpr Matrix<1,ZFA01::N_AX> Rg = ZFA01::Rg;			// [-] 減速比
	static constexpr Matrix<1,ZFA01::N_AX> Ktn = ZFA01::Kt;			// [Nm/A]  トルク定数
	static constexpr Matrix<1,ZFA01::N_AX> Jan = ZFA01::Ja;			// [kgm^2] 全慣性
	static constexpr Matrix<1,ZFA01::N_AX> gdis = {200, 200, 200};	// [rad/s] 外乱オブザーバの推定帯域
	static constexpr double wp = 30, zp = 1;						// [rad/s] 位置制御帯域, [-] 位置制御制動係数
	static constexpr size_t SrefN = 10000;							// [samples] S字軌道生成の窓幅
	static constexpr Matrix<1,ZFA01::N_AX> Kpv  = {2.0*zp*wp, 2.0*zp*wp, 2.0*zp*wp};	// [rad/s] 速度制御比例ゲイン(=速度制御帯域)
	static constexpr Matrix<1,ZFA01::N_AX> Kpp  = {wp/2.0/zp, wp/2.0/zp, wp/2.0/zp};	// [-] 位置制御比例ゲイン
	static constexpr Matrix<3,3> I = Matrix<3,3>::eye();			// 単位行列
	static constexpr double Pgra = 100e-3;							// [m] 研磨送り振幅指令
	static constexpr Matrix<1,3> pGsta = {325e-3, 2.00e-3, 118e-3};	// [m] 研磨開始位置
	
	// 制御用変数宣言
	static double Cfset = 0.05;	// (オンライン設定用) [1/kg] 力制御ゲイン
	static double Fdmp = 10;	// (オンライン設定用) [Hz] HPFダンピングのカットオフ周波数
	static double delta = 10;	// (オンライン設定用) [-]  HPFダンピングの粘性係数
	static double Fgra = -5;	// (オンライン設定用) [N] 研磨力指令
	static Matrix<1,ZFA01::N_AX> Cf = {Cfset, Cfset, Cfset};	// [1/kg]  力制御器 比例ゲイン
	static Matrix<1,ZFA01::N_AX> thm, thl;		// [rad] 位置ベクトル モータ側, 負荷側
	static Matrix<1,ZFA01::N_AX> wm, wl;		// [rad/s] 速度ベクトル モータ側, 負荷側
	static Matrix<1,ZFA01::N_AX> taus;			// [Nm] ねじれトルクベクトル
	static Matrix<1,ZFA01::N_AX> tdis;			// [Nm] 推定外乱トルクベクトル
	static Matrix<1,ZFA01::N_AX> amref, alref;	// [rad/s^2] 加速度指令ベクトル モータ側, 負荷側
	static Matrix<1,ZFA01::N_AX> iqref;			// [A] q軸電流指令ベクトル
	static Matrix<1,3> p1, p2, p3;				// [m] 作業空間位置ベクトル[px,py,pz]^T 1S軸, 2L軸, 3U軸
	static Matrix<1,3> pG, vG;					// [m] 作業空間位置ベクトル[px,py,pz]^T, 速度ベクトル[vx,vy,vz]^T 先端G(グラインド点)
	static Matrix<1,3> pGptp;					// [m] 作業空間位置指定(Point-to-Point) 先端G(グラインド点)
	static Matrix<1,3> pGref, vGref, aGref;		// [m] 作業空間位置指令, [m/s] 速度指令, [m/s^2] 加速度指令 先端G(グラインド点)
	static Matrix<1,3> fsens;					// [N] 先端6軸力覚センサの検出推力(XYZ並進分)
	static Matrix<1,3> fG;						// [N] 作業空間推力ベクトル[fx,fy,fz]^T 先端G(グラインド点)
	static Matrix<1,3> fGref;					// [N] 作業空間推力指令 先端G(グラインド点)
	static Matrix<3,3> S = {1,0,0, 0,1,0, 0,0,1};	// ハイブリッド制御 選択行列(1 = 位置制御, 0 = 力制御)
	static Matrix<1,3> aGreff, aGrefp;				// [m/s^2] 力制御側の加速度指令, 位置制御側の加速度指令
	
	static ZFA01 Robot;								// ZFA01ロボット制御クラス
	static DisturbanceObsrv<DObType::FULL_0TH, ZFA01::N_AX> DOB(Ktn, Jan, gdis, Ts);// 外乱オブザーバ(3軸分)
	static SshapeGenerators<SrefN, 3> PosRefShaper;									// S字軌道生成器(3軸分)
	static HighPassFilter HPF(2.0*M_PI*Fdmp, Ts);		// ダンピング用HPF
	
	if(CmdFlag == CTRL_INIT){
		// 初期化モード (ここは制御開始時/再開時に1度だけ呼び出される(非リアルタイム空間なので重い処理もOK))
		Initializing = true;				// 初期化中ランプ点灯
		Screen.InitOnlineSetVar(Cfset, Fdmp, delta, Fgra);	// オンライン設定変数の初期値の設定
		Interface.ServoON();				// サーボON指令の送出
		Interface.SetForceFilterLow();		// 力覚センサのフィルタを低域に設定
		
		// 位置指令系 初期化処理
		Interface.GetPosition(PositionRes);	// [rad] 初期位置の取得
		thm.LoadArray(PositionRes);			// [rad] モータ側初期位置
		thl = thm % Rg;						// [rad] 負荷側位置(簡易的換算)
		Robot.GetWorkspacePosition(thl, p1, p2, p3, pG);// [m] 作業空間位置の計算
		PosRefShaper.SetInitialValue(pG);	// S字軌道生成器の初期値を初期位置に設定
		
		Initializing = false;				// 初期化中ランプ消灯
	}
	if(CmdFlag == CTRL_LOOP){
		// 周期モード (ここは制御周期 SAMPLING_TIME[0] 毎に呼び出される(リアルタイム空間なので処理は制御周期内に収めること))
		// リアルタイム制御ここから
		Interface.GetPositionAndVelocity(PositionRes, VelocityRes);	// [rad,rad/s] 位置と速度応答の取得
		Interface.GetTorque(TorqueRes);						// [Nm] ねじれトルク応答の取得
		Screen.GetOnlineSetVar(Cfset, Fdmp, delta, Fgra);	// オンライン設定変数の読み込み
		
		// ゲインパラメータのオンライン設定
		Cf.Set(Cfset, Cfset, Cfset);	// 力制御ゲイン
		HPF.SetCutFreq(2.0*M_PI*Fdmp);	// HPFダンピングのカットオフ周波数
		
		// センサ系データを縦ベクトルとして読み込んで整理
		thm.LoadArray(PositionRes);		// [rad] モータ側位置
		wm.LoadArray(VelocityRes);		// [rad/s] モータ側速度
		taus.LoadArray(TorqueRes);		// [Nm] ねじれトルク
		thl = thm % Rg;					// [rad] 負荷側位置(簡易的換算)
		wl = wm % Rg;					// [rad/s] 負荷側速度(簡易的換算)
		fsens.LoadArray(ForceRes);		// [N] 先端G座標系の並進推力
		
		// 関節空間から作業空間への変換
		Robot.GetWorkspacePosition(thl, p1, p2, p3, pG);// [m] 作業空間位置の取得
		Robot.GetWorkspaceForce(fsens, fG);				// [N] 作業空間推力の取得
		Robot.UpdateJacobi(thl);		// Jacobi行列の更新
		vG = Robot.Jaco(wl);			// [m/s] 作業空間速度の計算
		
		// 作業空間力指令の生成
		if(0 <= t && t < 7){
			fGref = Matrix<1,3>::zeros();	// [N] 0～7秒は，ゼロ推力指令
		}else{
			fGref.Set(0, 0, Fgra);			// [N] 7秒以降は，Z方向研磨力指令
		}
		
		// 作業空間位置指令の生成
		if(0 <= t && t < 2){
			pGptp = ZFA01::pGhome;			// [m] 0～2秒はホームポジションに移動
		}else if(2 <= t && t < 3){
			Interface.SetZero6axisForce();	// 2～3秒は力覚センサのゼロキャリブレーション
		}else if(3 <= t && t < 6){
			pGptp = pGsta; 					// [m] 3～6秒は研磨ポジションへ移動
		}else if(6 <= t && t < 8){
											// 6～8秒は位置指令生成はしない
		}else{
			pGptp[1] = pGsta[1] + PulseWave(0.25, 0, t, 8)*Pgra;// [m] 8秒以降はX方向研磨送り指令
		}
		pGref = PosRefShaper.GetShapedSignal(pGptp);			// [m] 位置指令(S字成形後)
		
		// 力制御/位置制御の切り替え
		if(0 <= t && t < 6){
			S.Set(1,0,0, 0,1,0, 0,0,1);	// 0～6秒は，全方向位置制御
		}else{
			S.Set(1,0,0, 0,1,0, 0,0,0); // 6秒以降は，Z方向力制御，X,Y方向位置制御
		}
		
		// 力/位置ハイブリッド制御
		aGreff = (fGref - fG) & Cf;					// [m/s^2] 力制御側 加速度指令演算(作業空間)
		aGreff[3] = aGreff[3] - delta*HPF.GetSignal(vG[3]);	// HPFダンピング(作業空間Z方向)
		vGref  = (pGref - pG) & Kpp;				// [m/s]   位置制御側 速度指令演算(作業空間)
		aGrefp = (vGref - vG) & Kpv;				// [m/s^2] 位置制御側 加速度指令演算(作業空間)
		aGref = (I - S)*aGreff + S*aGrefp;			// [m/s^2] 加速度指令(作業空間)
		
		// 加速度制御
		alref = Robot.Jacoinv(aGref);				// [rad/s^2] 加速度指令変換(関節空間，負荷側)
		amref = alref & Rg;							// [rad/s^2] 加速度指令換算(関節空間，モータ側)
		tdis = DOB.GetDistTorque(iqref, wm);		// [Nm] 推定外乱演算
		iqref = (amref & Jan % Ktn) + (tdis % Ktn);	// [A]  電流指令
		
		// 安全系と電流指令設定
		Robot.CheckSafety(thl, pG);		// 安全監視（範囲外になると緊急停止）
		Robot.CurrentLimiter(iqref);	// [A] q軸電流リミッタ
		iqref.StoreArray(CurrentRef);	// [A] 最終的な全軸の電流指令
		
		Interface.SetCurrent(CurrentRef);	// [A] 電流指令の出力
		Screen.SetVarIndicator(pG[1], pG[2], pG[3], fG[1], fG[2], fG[3], pGref[1], pGref[2], fGref[3], 0);	// 任意変数インジケータ(変数0, ..., 変数9)
		Graph.SetTime(Tact, t);											// [s] グラフ描画用の周期と時刻のセット
		Graph.SetVars(0, iqref[1], iqref[2], iqref[3], 0, 0, 0, 0, 0);	// グラフプロット0 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(1, thm[1], thm[2], thm[3], 0, 0, 0, 0, 0);		// グラフプロット1 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(2, wm[1], wm[2], wm[3], 0, 0, 0, 0, 0);			// グラフプロット2 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(3, taus[1], taus[2], taus[3], 0, 0, 0, 0, 0);		// グラフプロット3 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(4, vG[1], vG[2], vG[3], 0, 0, 0, 0, 0);			// グラフプロット4 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetVars(5, fG[1], fG[2], fG[3], 0, 0, 0, 0, 0);			// グラフプロット5 (グラフ番号, 変数0, ..., 変数7)
		Graph.SetWorkspace(p2, p3, pG);									// 作業空間グラフプロット
		Memory.SetData(Tact, t, pG[1], pG[2], pG[3], fG[1], fG[2], fG[3], pGref[1], pGref[2], fGref[3]);	// CSVデータ保存変数 (周期, A列, B列, ..., J列)
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
		Interface.Get6axisForce(ForceRes);	// [N,Nm] 6軸力覚応答の取得
		
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

