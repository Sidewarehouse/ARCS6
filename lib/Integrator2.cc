//! @file Integrator2.cc
//! @brief 2次積分器クラス
//!
//! 2次積分器 G(s)=1/s^2 (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "Integrator2.hh"

using namespace ARCS;

Integrator2::Integrator2(double SmplTime)
	// コンストラクタ SmplTime；[s] 制御周期
	: Ts(SmplTime),		// [s] 制御周期の格納
	  uZ1(0), uZ2(0), yZ1(0), yZ2(0)
{
	
}

Integrator2::~Integrator2(){
	// デストラクタ
}

double Integrator2::GetSignal(double u){
	// 出力信号の取得 u；入力信号
	double y;
	
	y=(Ts*Ts)/4.0* ( u + 2.0*uZ1 + uZ2) + 2.0*yZ1 - yZ2;

	uZ2=uZ1;
	uZ1=u;
	yZ2=yZ1;
	yZ1=y;

	return y;
}

void Integrator2::SetSmplTime(double SmplTime){
	// 制御周期の再設定 SmplTime；[s] 制御周期
	Ts=SmplTime;	// [s] 制御周期の再設定
}

void Integrator2::ClearStateVars(void){
	// すべての状態変数のリセット
	uZ1=0;	// 状態変数1のゼロクリア
	uZ2=0;	// 状態変数2のゼロクリア
	yZ1=0;	// 状態変数3のゼロクリア
	yZ2=0;	// 状態変数4のゼロクリア
}


