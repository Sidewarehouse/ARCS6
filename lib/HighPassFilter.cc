//! @file HighPassFilter.cc
//! @brief 高域通過濾波器クラス
//!
//! 1次高域通過濾波器 G(s)=s/(s+g) (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "HighPassFilter.hh"

using namespace ARCS;

HighPassFilter::HighPassFilter(double CutFreq, double SmplTime)
	// コンストラクタ CutFreq；[rad/s] 遮断周波数，SmplTime；[s] 制御周期
	: Ts(SmplTime),	// [s]		制御周期の格納
	  g(CutFreq),	// [rad/s]	遮断周波数
	  uZ1(0), yZ1(0)
{
	
}

HighPassFilter::~HighPassFilter(){
	// デストラクタ
}

double HighPassFilter::GetSignal(double u){
	// 出力信号の取得 u；入力信号
	double y;
	
	y=1.0/(g*Ts+2.0)*( 2.0*(u-uZ1) - yZ1*(g*Ts-2.0) );
	
	uZ1=u;
	yZ1=y;
	
	return y;
}

void HighPassFilter::SetCutFreq(double CutFreq){
	// 遮断周波数の再設定 CutFreq；[rad/s] 遮断周波数
	g=CutFreq;		// [rad/s]	遮断周波数の再設定
}

void HighPassFilter::SetSmplTime(double SmplTime){
	// 制御周期の再設定 SmplTime；[s] 制御周期
	Ts=SmplTime;	// [s] 制御周期の再設定
}

void HighPassFilter::ClearStateVars(void){
	// すべての状態変数のリセット
	uZ1=0;	// 状態変数1のゼロクリア
	yZ1=0;	// 状態変数2のゼロクリア
}



