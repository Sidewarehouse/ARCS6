//! @file MotorSimulator.cc
//! @brief モータシミュレータ
//!
//! モータを模擬する
//!
//! @date 2020/04/14
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include <tuple>
#include "MotorSimulator.hh"
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

//! @brief 空コンストラクタ
MotorSimulator::MotorSimulator(void)
	: Kt(0), Jm(0), Dm(0), Ts(0), iq(0), taul(0),
	  A(), B(), u(), y(), PlantSys()
{
	PassedLog();
}

//! @brief コンストラクタ
//! @param[in]	TrqConst	[Nm/A] トルク定数
//! @param[in]	MotorInert	[kgm^2]モータ慣性
//! @param[in]	MotorFric	[Nm/(rad/s)]モータ粘性
//! @param[in]	SmplTime	[s] サンプリング時間
MotorSimulator::MotorSimulator(const double TrqConst, const double MotorInert, const double MotorFric, const double SmplTime)
	: Kt(TrqConst), Jm(MotorInert), Dm(MotorFric), Ts(SmplTime), iq(0), taul(0),
	  A(), B(), u(), y(), PlantSys()
{
	SetStateSpaceModel();	// 状態空間モデルにセット
	PassedLog();
}

//! @brief ムーブコンストラクタ
MotorSimulator::MotorSimulator(MotorSimulator&& r)
	: Kt(r.Kt), Jm(r.Jm), Dm(r.Dm), Ts(r.Ts), iq(r.iq), taul(r.taul),
	  A(r.A), B(r.B), u(r.u), y(r.y), PlantSys(std::move(r.PlantSys))
{
	
}

//! @brief デストラクタ
MotorSimulator::~MotorSimulator(){
	PassedLog();
}

//! @brief パラメータを設定する関数
//! @param[in]	TrqConst	[Nm/A] トルク定数
//! @param[in]	MotorInert	[kgm^2]モータ慣性
//! @param[in]	MotorFric	[Nm/(rad/s)]モータ粘性
//! @param[in]	SmplTime	[s] サンプリング時間
void MotorSimulator::SetParameters(const double TrqConst, const double MotorInert, const double MotorFric, const double SmplTime){
	Kt = TrqConst;
	Jm = MotorInert;
	Dm = MotorFric;
	Ts = SmplTime;
	SetStateSpaceModel();	// 状態空間モデルにセット
}

//! @brief 予めモータ電流と負荷トルクを設定する関数
//! @param[in]	current	[A] 電流
//! @param[in]	loadtorque	[Nm] 負荷トルク
void MotorSimulator::SetCurrentAndLoadTorque(const double current, const double loadtorque){
	iq = current;
	taul = loadtorque;
}

//! @brief モータ速度と位置を返す関数(普通版)
//! @param[in]	current	[A] 電流
//! @param[in]	loadtorque	[Nm] 負荷トルク
//! @param[out]	velocity	モータ速度 [rad/s]
//! @param[out]	position	位置 [rad]
void MotorSimulator::GetResponses(const double current, const double loadtorque, double& velocity, double& position){
	// 入力ベクトルの設定
	u.Set(
		current,
		loadtorque
	);
	
	// 応答計算
	PlantSys.GetResponses(u, y);
	
	// 状態ベクトルから抽出
	velocity = y[1];
	position = y[2];
}

//! @brief モータ速度と位置を返す関数(タプル版)
//! @param[in]	current	電流 [A]
//! @param[in]	loadtorque	負荷トルク [Nm]
//! @return	モータ速度 [rad/s], 位置 [rad]
std::tuple<double,double> MotorSimulator::GetResponses(const double current, const double loadtorque){
	double velocity, position;
	GetResponses(current, loadtorque, velocity, position);	// 応答計算
	
	// 速度と位置をタプルで返す
	return std::forward_as_tuple(velocity, position);
}

//! @brief 予めセットされたモータ電流と負荷トルクから速度と位置を返す関数(普通版)
//! @param[out]	velocity	速度 [rad/s]
//! @param[out] position	位置 [rad]
void MotorSimulator::GetVelocityAndPosition(double& velocity, double& position){
	GetResponses(iq, taul, velocity, position);
}

//! @brief 予めセットされたモータ電流と負荷トルクから速度と位置を返す関数(タプル版)
//! @return モータ速度 [rad/s]，位置 [rad]
std::tuple<double, double> MotorSimulator::GetVelocityAndPosition(void){
	return GetResponses((const double)iq, (const double)taul);
}

//! @brief モータ慣性を設定する関数
//! @param[in]	inertia	モータ慣性 [kgm^2]
void MotorSimulator::SetMotorInertia(const double inertia){
	Jm = inertia;
	SetStateSpaceModel();	// 状態空間モデルにセット
}

//! @brief シミュレータをリセットする関数
void MotorSimulator::Reset(void){
	PlantSys.ClearStateVector();
}

//! @brief 状態空間モデルを設定する関数
void MotorSimulator::SetStateSpaceModel(void){
	// 連続系A行列の設定
	A.Set(
		-Dm/Jm,  0,
		     1,  0
	);
	
	// 連続系B行列の設定
	B.Set(
		Kt/Jm, -1.0/Jm,
		    0,       0
	);
	
	auto C = Matrix<2,2>::eye();		// C行列の設定
	PlantSys.SetContinuous(A, B, C, Ts);// 状態空間モデルの設定
}

