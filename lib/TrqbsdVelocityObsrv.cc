//! @file TrqbsdVelocityObsrv.cc
//! @brief トルクセンサベース速度オブザーバ
//! @date 2018/12/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2018 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "TrqbsdVelocityObsrv.hh"
#include <math.h>

using namespace ARCS;

//! @brief コンストラクタ
TrqbsdVelocityObsrv::TrqbsdVelocityObsrv(
		double TrqConst, double Inertia, double Viscosity, double GearRatio, double Bandwidth, double SmplTime
	)
	: Ktn(TrqConst), Jmn(Inertia), Dmn(Viscosity), Rg(GearRatio), w(Bandwidth), Ts(SmplTime),
	  e(exp(-Dmn/Jmn*Ts)),					// e = exp(-Dmn/Jmn*Ts)
	  k1(Dmn/Jmn*pow(exp(-w*Ts)-e,2)/(1-e)),// オブザーバゲイン1の計算
	  k2(1 + e - 2*exp(-w*Ts)),				// オブザーバゲイン2の計算
	  k(),				// オブザーバのゲインベクトル
	  Ad(), Bd(), c(),	// 離散状態方程式＆出力方程式の係数行列
	  u(), xh(), xh_next(), y()	// 信号ベクトル
{
	// プラントの離散系状態方程式のA行列
	Ad.Set(
		e            , 0,
		Jmn/Dmn*(1-e), 1
	);
	
	// プラントの離散系状態方程式のB行列
	Bd.Set(
		Ktn/Dmn*(1-e)              , -1/(Dmn*Rg)*(1-e)             ,
		Ktn/Dmn*(Ts-Jmn/Dmn*(1-e)) , -1/(Dmn*Rg)*(Ts-Jmn/Dmn*(1-e))
	);
	
	// プラントの出力方程式
	c.Set(0, 1);
	
	// オブザーバゲインベクトル
	k.Set(
		k1,
		k2
	);
}

//! @brief デストラクタ
TrqbsdVelocityObsrv::~TrqbsdVelocityObsrv(){
	
}

//! @brief 速度推定値を計算する関数
double TrqbsdVelocityObsrv::GetVelocity(double Current, double TorsionTorque, double Position){
	double Velocity = 0;
	
	// プラント入力ベクトルの設定
	u.Set(
		Current,
		TorsionTorque
	);
	
	// プラント出力ベクトルの設定
	y.Set(Position);
	
	// 離散状態方程式の計算
	xh_next = Ad*xh + Bd*u + k*(y - c*xh);
	
	// 状態推定ベクトルから速度を抽出
	Velocity = xh.GetElement(1,1);
	
	// 状態推定ベクトルの更新
	xh = xh_next;
	
	return Velocity;
}



