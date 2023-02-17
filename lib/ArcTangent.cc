//! @file ArcTangent.cc
//! @brief アークタンジェント(いわゆるatan2を拡張した「atan3」のようなもの)
//!
//! 普通のatan2は±πの範囲しか使えないが，これは全領域で使える。つまり多回転もオーケー。
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include <cmath>
#include "ArcTangent.hh"

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
ArcTangent::ArcTangent(void)
	: RotationCounter(0),
	  theta_z1(0)
{
	PassedLog();
}

//! @brief コンストラクタ(初期角度を与える場合)
//! @param[in] theta_ini 初期角度 [rad]
ArcTangent::ArcTangent(const double theta_ini)
	: RotationCounter(0),
	  theta_z1(theta_ini)
{
	PassedLog();
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	右辺値
ArcTangent::ArcTangent(ArcTangent&& r)
	: RotationCounter(r.RotationCounter),
	  theta_z1(r.theta_z1)
{
	
}

//! @brief デストラクタ
ArcTangent::~ArcTangent(){
	PassedLog();
}

//! @brief アークタンジェントを計算して角度を返す関数「atan3」
//! @param[in]	Num	atan3の分子
//! @param[in]	Den	atan3の分母
//! @return tan(Num/Den)の角度 [rad]
double ArcTangent::GetAngle(const double Num, const double Den){
	const double thbuff = atan2(Num, Den);
	if((thbuff < -M_PI_2) && (M_PI_2 < theta_z1)){
		// 正から負に移動したら，+1回転したとみなす
		++RotationCounter;
	}else if((theta_z1 < -M_PI_2) && (M_PI_2 < thbuff)){
		// 負から正に移動したら，-1回転したとみなす
		--RotationCounter;
	}
	theta_z1 = thbuff;	// 次回用に取っておく
	return thbuff + (double)RotationCounter*2.0*M_PI;	// atan2の結果に回転回数分のradを加算して返す
}

//! @brief 角度を指定した値で初期化する関数
//! @param[in]	theta_ini	初期角度 [rad]
void ArcTangent::InitAngle(double theta_ini){
	theta_z1 = theta_ini;
}

//! @brief 回転回数カウンタを設定する関数
//! @param[in]	rot_init	回転回数 [-]
void ArcTangent::InitRotation(const long rot_init){
	RotationCounter = rot_init;
}
