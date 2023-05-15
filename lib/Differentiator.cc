//! @file Differentiator.cc
//! @brief  擬似微分器クラス
//!
//! 擬似微分器 G(s)=(s*gpd)/(s+gpd) (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "Differentiator.hh"

using namespace ARCS;

Differentiator::Differentiator(double Bandwidth, double SmplTime)
	// コンストラクタ Bandwidth；[rad/s] 帯域，SmplTime；[s] 制御周期
	: Ts(SmplTime),		// [s]		制御周期の格納
	  gpd(Bandwidth),	// [rad/s]	擬似微分の帯域の格納
	  uZ1(0), yZ1(0)
{
	
}

Differentiator::~Differentiator(){
	// デストラクタ
}

double Differentiator::GetSignal(double u){
	// 出力信号の取得 u；入力信号
	double y;
	
	y = ( 2.0*gpd*(u-uZ1) + (2.0-Ts*gpd)*yZ1 )/(2.0+Ts*gpd);
	
	uZ1=u;
	yZ1=y;
	
	return y;
}

void Differentiator::SetBandwidth(double Bandwidth){
	// 擬似微分の帯域の再設定 Bandwidth；[rad/s] 帯域
	gpd=Bandwidth;
}

void Differentiator::SetSmplTime(double SmplTime){
	// 制御周期の再設定 SmplTime；[s] 制御周期
	Ts=SmplTime;	// [s] 制御周期の再設定
}

void Differentiator::ClearStateVars(void){
	// すべての状態変数のリセット
	uZ1=0;	// 状態変数1のゼロクリア
	yZ1=0;	// 状態変数2のゼロクリア
}


