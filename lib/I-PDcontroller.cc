//! @file I-PDcontroller.cc
//! @brief I-PD制御器
//! @date 2020/03/05
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include <tuple>
#include "I-PDcontroller.hh"
#include "Discret.hh"

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
I_PDcontroller::I_PDcontroller(double Igain, double Pgain, double Dgain, double Bandwidth, double SmplTime)
	: Ki(Igain), Kp(Pgain), Kd(Dgain), g(Bandwidth), Ts(SmplTime),
	  Ac(), Bc(), Ad(), Bd(), c(), d(), u(), x(), x_next(), y()
{
	// 連続系状態方程式のシステム行列
	Ac.Set(
		0,	 0,
		0,	-g
	);
	Bc.Set(
		1,	-1,
		0,	 g
	);
	c.Set(
		Ki,	g*Kd
	);
	d.Set(
		0,	-Kp - g*Kd
	);
	
	// 離散化
	std::tie(Ad,Bd) = Discret::GetDiscSystem(Ac,Bc,Ts);
	
	PassedLog();
}

//! @brief デストラクタ
I_PDcontroller::~I_PDcontroller(){
	PassedLog();
}

//! @brief 制御器の出力信号を計算して出力する関数
//! @param[in]	ref	指令値
//! @param[in]	res	応答値
//! @return	出力信号
double I_PDcontroller::GetSignal(double ref, double res){
	
	// 入力ベクトルの設定
	u.Set(
		ref,
		res
	);
	
	// 離散系状態方程式の計算
	x_next = Ad*x + Bd*u;
	y = c*x + d*u;
	
	// 状態ベクトルの更新
	x = x_next;

	// 出力ベクトルから抽出
	return y.GetElement(1,1);
}

//! @brief 積分ゲインを再設定する関数
void I_PDcontroller::SetIgain(double Igain){
	Ki = Igain;
}

//! @brief 先行比例ゲインを再設定する関数
void I_PDcontroller::SetPgain(double Pgain){
	Kp = Pgain;
}

//! @brief 先行微分ゲインを再設定する関数
void I_PDcontroller::SetDgain(double Dgain){
	Kd = Dgain;
}

//! @brief 擬似微分の帯域を再設定する関数
void I_PDcontroller::SetBandwidth(double Bandwidth){
	g = Bandwidth;
}

