//! @file PIDcontroller.cc
//! @brief PID制御器クラス
//!
//! PID制御器 G(s) = Kp + Ki/s + Kd*s*gdis/(s+gdis) (双一次変換)
//!
//! @date 2019/07/29
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "PIDcontroller.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

using namespace ARCS;

//! @brief コンストラクタ
//! @param[in]	Pgain	比例ゲイン
//! @param[in]	Igain	積分ゲイン
//! @param[in]	Dgain	微分ゲイン
//! @param[in]	SmplTime	[s] 制御周期
PIDcontroller::PIDcontroller(double Pgain, double Igain, double Dgain, double Bandwidth, double SmplTime)
	: Ts(SmplTime),		// [s]		制御周期の格納
	  Kp(Pgain),		// 			比例ゲインの格納
	  Ki(Igain),		// 			積分ゲインの格納
	  Kd(Dgain),		// 			微分ゲインの格納
	  gpd(Bandwidth),	// [rad/s]	擬似微分の帯域の格納
	  uZ1(0), uZ2(0), yZ1(0), yZ2(0)
{
	PassedLog();
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	右辺値
PIDcontroller::PIDcontroller(PIDcontroller&& r)
	: Ts(r.Ts),		// [s]		制御周期の格納
	  Kp(r.Kp),		// 			比例ゲインの格納
	  Ki(r.Ki),		// 			積分ゲインの格納
	  Kd(r.Kd),		// 			微分ゲインの格納
	  gpd(r.gpd),	// [rad/s]	擬似微分の帯域の格納
	  uZ1(r.uZ1), uZ2(r.uZ2), yZ1(r.yZ1), yZ2(r.yZ2)
{
	
}

//! @brief デストラクタ
PIDcontroller::~PIDcontroller(){
	PassedLog();
}

//! @brief 出力信号の取得
//! @param[in]	u	PID制御器入力
//! @return	PID制御器出力
double PIDcontroller::GetSignal(double u){
	double y = ( 2.0*Ts*(Ki+Kp*gpd)*(u-uZ2) + Ki*gpd*Ts*Ts*(u+2.0*uZ1+uZ2) + 4.0*(Kd*gpd+Kp)*(u-2.0*uZ1+uZ2) - (4.0-2.0*gpd*Ts)*yZ2 + 8.0*yZ1 )/(4.0+2.0*gpd*Ts);
	
	uZ2 = uZ1;
	uZ1 = u;
	yZ2 = yZ1;
	yZ1 = y;
	
	return y;
}

//! @brief 比例ゲインの再設定
//! @param[in]	Pgain	比例ゲイン
void PIDcontroller::SetPgain(double Pgain){
	Kp = Pgain;
}

//! @brief 積分ゲインの再設定
//! @param[in]	Igain	積分ゲイン
void PIDcontroller::SetIgain(double Igain){
	Ki = Igain;
}

//! @brief 微分ゲインの再設定
//! @param[in]	Dgain	微分ゲイン
void PIDcontroller::SetDgain(double Dgain){
	Kd = Dgain;
}

//! @brief 擬似微分の帯域の再設定
//! @param[in]	Bandwidth	[rad/s] 帯域
void PIDcontroller::SetBandwidth(double Bandwidth){
	gpd = Bandwidth;
}

//! @brief 制御周期の再設定
//! @param[in]	SmplTime	[s] 制御周期
void PIDcontroller::SetSmplTime(double SmplTime){
	Ts = SmplTime;	// [s] 制御周期の再設定
}

//! @brief すべての状態変数のリセット
void PIDcontroller::ClearStateVars(void){
	uZ1 = 0;	// 状態変数1のゼロクリア
	uZ2 = 0;	// 状態変数2のゼロクリア
	yZ1 = 0;	// 状態変数3のゼロクリア
	yZ2 = 0;	// 状態変数4のゼロクリア
}

