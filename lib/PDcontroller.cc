//! @file PDcontroller.cc
//! @brief PD制御器クラス
//!
//! PD制御器 G(s) = Kp + Kd*s*gdis/(s+gdis) (双一次変換)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "PDcontroller.hh"

using namespace ARCS;

PDcontroller::PDcontroller(double Pgain, double Dgain, double Bandwidth, double SmplTime)
	// コンストラクタ
	// Pgain；比例ゲイン，Dgain；微分ゲイン，SmplTime；[s] 制御周期
	: Ts(SmplTime),		// [s]		制御周期の格納
	  Kp(Pgain),		// 			比例ゲインの格納
	  Kd(Dgain),		// 			微分ゲインの格納
	  gpd(Bandwidth),	// [rad/s]	擬似微分の帯域の格納
	  uZ1(0), yZ1(0)
{
	
}

PDcontroller::~PDcontroller(){
	// デストラクタ
}

double PDcontroller::GetSignal(double u){
	// 出力信号の取得 u；入力信号
	double y;
	
	y = (  (2.0-gpd*Ts)*yZ1 + Kp*( 2.0*(u-uZ1) + gpd*Ts*(u+uZ1) ) + 2.0*Kd*gpd*(u-uZ1) )/(2.0+gpd*Ts);
	
	uZ1=u;
	yZ1=y;
	
	return y;
}

void PDcontroller::SetPgain(double Pgain){
	// 比例ゲインの再設定 Pgain；比例ゲイン
	Kp=Pgain;
}

void PDcontroller::SetDgain(double Dgain){
	// 微分ゲインの再設定 Dgain；微分ゲイン
	Kd=Dgain;
}

void PDcontroller::SetBandwidth(double Bandwidth){
	// 擬似微分の帯域の再設定 Bandwidth；[rad/s] 帯域
	gpd=Bandwidth;
}

void PDcontroller::SetSmplTime(double SmplTime){
	// 制御周期の再設定 SmplTime；[s] 制御周期
	Ts=SmplTime;	// [s] 制御周期の再設定
}

void PDcontroller::ClearStateVars(void){
	// すべての状態変数のリセット
	uZ1=0;	// 状態変数1のゼロクリア
	yZ1=0;	// 状態変数2のゼロクリア
}


