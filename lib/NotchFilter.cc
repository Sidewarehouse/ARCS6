//! @file NotchFilter.cc
//! @brief ノッチフィルタクラス
//!
//! ノッチフィルタ G(s)=( s^2 + w^2 )/( s^2 + w/Q*s + w^2) (双一次変換)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "NotchFilter.hh"

using namespace ARCS;

NotchFilter::NotchFilter(double AngFreq, double Sharp, double SmplTime)
	// コンストラクタ AngFreq；[rad/s] 遮断中心周波数，Sharp；鋭さ，SmplTime；[s] 制御周期
	: Ts(SmplTime),	// [s]		制御周期の格納
	  w(AngFreq),	// [rad/s]	遮断中心周波数の格納
	  Q(Sharp),		//			鋭さの格納
	  uZ1(0), uZ2(0), yZ1(0), yZ2(0)
{
	
}

NotchFilter::~NotchFilter(){
	// デストラクタ
}

double NotchFilter::GetSignal(double u){
	// 出力信号の取得 u；入力
	double y;
	
	y=1.0/(4.0+2.0*Ts*w/Q+w*w*Ts*Ts)*( (u+uZ2)*(4.0+w*w*Ts*Ts) + (uZ1-yZ1)*(2.0*w*w*Ts*Ts-8.0) - yZ2*(4.0-2.0*Ts*w/Q+w*w*Ts*Ts) );

	uZ2=uZ1;
	uZ1=u;
	yZ2=yZ1;
	yZ1=y;
	
	return y;
}

void NotchFilter::SetCuttoff(double AngFreq){
	// 遮断中心周波数の再設定 AngFreq；[rad/s] 帯域
	w=AngFreq;
}

void NotchFilter::SetSharpness(double Sharp){
	// 鋭さの再設定 Sharp；鋭さ
	Q=Sharp;
}

void NotchFilter::SetSmplTime(double SmplTime){
	// 制御周期の再設定 SmplTime；[s] 制御周期
	Ts=SmplTime;	// [s] 制御周期の再設定
}

void NotchFilter::ClearStateVars(void){
	// すべての状態変数のリセット
	uZ1=0;	// 状態変数1のゼロクリア
	uZ2=0;	// 状態変数2のゼロクリア
	yZ1=0;	// 状態変数3のゼロクリア
	yZ2=0;	// 状態変数4のゼロクリア
}


