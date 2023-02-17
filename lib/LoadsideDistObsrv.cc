//! @file LoadsideDistObsrv.cc
//! @brief 同一次元負荷側外乱オブザーバ
//!
//! モータ側速度と捻れトルクから負荷側外乱トルクを推定します。
//!
//! @date 2018/01/17
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <tuple>
#include "LoadsideDistObsrv.hh"
#include "Discret.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

//! @brief コンストラクタ(オブザーバの極を3重根にする場合)
//! @param [in] Stiffness 2慣性間の剛性 [Nm/rad]
//! @param [in] LoadInertia 負荷側慣性 [kgm^2]
//! @param [in] GearRatio 減速比 [-]
//! @param [in] Bandwidth 推定帯域 [rad/s]
//! @param [in] SmplTime 制御周期 [s]
LoadsideDistObsrv::LoadsideDistObsrv(double Stiffness, double LoadInertia, double GearRatio, double Bandwidth, double SmplTime)
	: Ksn(Stiffness), Jln(LoadInertia), Rg(GearRatio), l1(-Bandwidth), l2(-Bandwidth), l3(-Bandwidth), Ts(SmplTime),
	  Ac(), Bc(), Ad(), Bd(), c(), u(), x(), x_next(), y()
{
	// 連続系A行列の設定
	Ac.Set(
					  l1 + l2 + l3, -Ksn,        0,
		(l3*(l1 + l2) + l1*l2)/Ksn,    0, -1.0/Jln,
				(Jln*l1*l2*l3)/Ksn,    0,        0
	);
	// 連続系B行列の設定
	Bc.Set(
		Ksn/Rg,                       - l1 - l2 - l3,
			 0, 1.0/Jln - (l3*(l1 + l2) + l1*l2)/Ksn,
			 0,                  -(Jln*l1*l2*l3)/Ksn
	);
	// C行列の設定
	c .Set(
		0, 0, 1
	);
	// 離散化
	std::tie(Ad,Bd) = Discret::GetDiscSystem(Ac,Bc,Ts);
	
	PassedLog();
}

//! @brief デストラクタ
LoadsideDistObsrv::~LoadsideDistObsrv(){
	
}

//! @brief 負荷側トルクを推定する関数
//! @param [in] MotorSpeed モータ側速度 [rad/s]
//! @param [in] TorsionTorque ねじれトルク [Nm]
double LoadsideDistObsrv::GetLoadsideTorque(double MotorSpeed, double TorsionTorque){
	// 入力ベクトルの設定
	u.Set(
		MotorSpeed,
		TorsionTorque
	);
	// 離散系状態方程式の計算
	x_next = Ad*x + Bd*u;
	y = c*x;
	// 状態ベクトルの更新
	x = x_next;
	// 出力ベクトルから抽出して返す
	return y.GetElement(1,1);
}


