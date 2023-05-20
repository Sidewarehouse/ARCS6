//! @file PhaseLeadLag.cc
//! @brief 位相進み/遅れ補償器クラス
//!
//! 位相補償器 G(s)=(1+s*a/g)/(1+s/g) (双一次変換)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "PhaseLeadLag.hh"

using namespace ARCS;

PhaseLeadLag::PhaseLeadLag(double AngFreq, double Phase, double SmplTime)
	// コンストラクタ
	// AngFreq；[rad/s] ゲイン特性の折れ点周波数，a < 1 で位相遅れ，1 < a で位相進み，SmplTime；[s] 制御周期
	: Ts(SmplTime),		// [s]		制御周期の格納
	  g(AngFreq),			// [rad/s]	折れ点周波数の格納
	  a(Phase),			//			a < 1 で位相遅れ，1 < a で位相進み
	  uZ1(0), yZ1(0)
{
	
}

PhaseLeadLag::~PhaseLeadLag(){
	// デストラクタ
}

double PhaseLeadLag::GetSignal(double u){
	// 出力信号の取得 u；入力信号
	double y;
	
	y=1.0/(Ts+2.0/g)*( u*(Ts+2.0*a/g) + uZ1*(Ts-2.0*a/g) - yZ1*(Ts-2.0/g) );
	
	uZ1=u;
	yZ1=y;
	
	return y;
}

void PhaseLeadLag::SetSmplTime(double SmplTime){
	// 制御周期の再設定 SmplTime；[s] 制御周期
	Ts=SmplTime;	// [s] 制御周期の再設定
}

void PhaseLeadLag::SetFreq(double AngFreq){
	// ゲイン特性の折れ点周波数の再設定 AngFreq；[rad/s] 折れ点周波数
	g=AngFreq;	// [rad/s] 折れ点周波数の再設定
}

void PhaseLeadLag::SetPhase(double Phase){
	// a < 1 で位相遅れ，1 < a で位相進み
	a=Phase;
}

void PhaseLeadLag::ClearStateVars(void){
	// すべての状態変数のリセット
	uZ1=0;	// 状態変数1のゼロクリア
	yZ1=0;	// 状態変数2のゼロクリア
}


