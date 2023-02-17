//! @file Differentiator2.cc
//! @brief 2次擬似微分器クラス
//!
//! 2次擬似微分器 G(s)=s^2*w^2/(s^2 + w/Q*s + w^2) (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "Differentiator2.hh"

using namespace ARCS;

Differentiator2::Differentiator2(double Bandwidth, double Sharp, double SmplTime)
	// コンストラクタ Bandwidth；[rad/s] 帯域，Sharp；鋭さ，SmplTime；[s] 制御周期
	: Ts(SmplTime),		// [s]		制御周期の格納
	  w(Bandwidth),		// [rad/s]	擬似微分の帯域の格納
	  Q(Sharp),			//			鋭さの格納
	  uZ1(0), uZ2(0), uZ3(0), yZ1(0), yZ2(0), yZ3(0)
{
	
}

Differentiator2::~Differentiator2(){
	// デストラクタ
}

double Differentiator2::GetSignal(double u){
	// 出力信号の取得 u；入力
	double y;
	
	y = 4.0*Q*Ts*w*w*(uZ3-uZ2-uZ1+u)
	   -4.0*Q*Ts*(yZ3-yZ2-yZ1)
	   -2.0*w*Ts*Ts*(-yZ3-yZ2+yZ1)
	   -Q*Ts*Ts*Ts*w*w*(yZ3+3.0*yZ2+3.0*yZ1);
	y = y * 1.0/(4.0*Q*Ts + 2.0*w*Ts*Ts + Q*Ts*Ts*Ts*w*w);
	
	uZ3=uZ2;
	uZ2=uZ1;
	uZ1=u;
	yZ3=yZ2;
	yZ2=yZ1;
	yZ1=y;
	
	return y;
}

void Differentiator2::SetBandwidth(double Bandwidth){
	// 擬似微分の帯域の再設定 Bandwidth；[rad/s] 帯域
	w=Bandwidth;
}

void Differentiator2::SetSharpness(double Sharp){
	// 鋭さの再設定 Sharp；鋭さ
	Q=Sharp;
}

void Differentiator2::SetSmplTime(double SmplTime){
	// 制御周期の再設定 SmplTime；[s] 制御周期
	Ts=SmplTime;	// [s] 制御周期の再設定
}

void Differentiator2::ClearStateVars(void){
	// すべての状態変数のリセット
	uZ1=0;	// 状態変数1のゼロクリア
	uZ2=0;	// 状態変数2のゼロクリア
	uZ3=0;	// 状態変数3のゼロクリア
	yZ1=0;	// 状態変数4のゼロクリア
	yZ2=0;	// 状態変数5のゼロクリア
	yZ3=0;	// 状態変数6のゼロクリア
}


