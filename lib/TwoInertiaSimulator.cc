//! @file TwoInertiaSimulator.hh
//! @brief 2慣性共振系シミュレータV4
//!
//! 2慣性共振系を模擬して，負荷側速度，ねじれ角，モータ速度/位置を計算して出力する。
//!
//! @date 2023/10/27
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <cassert>
#include <tuple>
#include "TwoInertiaParamDef.hh"
#include "TwoInertiaSimulator.hh"

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
TwoInertiaSimulator::TwoInertiaSimulator(void)
	: Jl(), Dl(), Ds(), Ks(),
	  Jm(), Dm(), Rg(), Kt(), Ts(), thl(), thm(),
	  PlantSys(), MotorInteg(), LoadInteg()
{
	PassedLog();
}

//! @brief コンストラクタ(Ds無し版)
//! @param[in]	Params	2慣性共振系パラメータ構造体
//! @param[in]	SmplTime	サンプリング周期 [s]
TwoInertiaSimulator::TwoInertiaSimulator(const struct TwoInertiaParamDef& Params, const double SmplTime)
	: Jl(Params.Jl), Dl(Params.Dl), Ds(0), Ks(Params.Ks),
	  Jm(Params.Jm), Dm(Params.Dm), Rg(Params.Rg), Kt(Params.Kt), Ts(SmplTime), thl(0), thm(0),
	  PlantSys(), MotorInteg(SmplTime), LoadInteg(SmplTime)
{
	SetStateSpaceModel();	// 状態空間モデルにセット
	PassedLog();
}

//! @brief コンストラクタ(Ds有り版)
//! @param[in]	Params		2慣性共振系パラメータ構造体
//! @param[in]	SmplTime	サンプリング周期 [s]
TwoInertiaSimulator::TwoInertiaSimulator(const struct TwoInertiaParams& Params, const double SmplTime)
	: Jl(Params.Jl), Dl(Params.Dl), Ds(Params.Ds), Ks(Params.Ks),
	  Jm(Params.Jm), Dm(Params.Dm), Rg(Params.Rg), Kt(Params.Kt), Ts(SmplTime), thl(0), thm(0),
	  PlantSys(), MotorInteg(SmplTime), LoadInteg(SmplTime)
{
	SetStateSpaceModel();	// 状態空間モデルにセット
	PassedLog();
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	右辺値
TwoInertiaSimulator::TwoInertiaSimulator(TwoInertiaSimulator&& r)
	: Jl(r.Jl), Dl(r.Dl), Ds(r.Ds), Ks(r.Ks),
	  Jm(r.Jm), Dm(r.Dm), Rg(r.Rg), Kt(r.Kt), Ts(r.Ts), thl(r.thl), thm(r.thm),
	  PlantSys(std::move(r.PlantSys)), MotorInteg(std::move(r.MotorInteg)), LoadInteg(std::move(r.LoadInteg))
{
	
}

//! @brief デストラクタ
TwoInertiaSimulator::~TwoInertiaSimulator(){
	PassedLog();
}

//! @brief 2慣性共振系のパラメータを設定する関数(Ds無し版)
//! @param[in]	Params		2慣性共振系パラメータ構造体
//! @param[in]	SmplTime	サンプリング周期 [s]
void TwoInertiaSimulator::SetParameters(const struct TwoInertiaParamDef& Params, const double SmplTime){
	Jl = Params.Jl;
	Dl = Params.Dl;
	Ds = 0;
	Ks = Params.Ks;
	Jm = Params.Jm;
	Dm = Params.Dm;
	Rg = Params.Rg;
	Kt = Params.Kt;
	Ts = SmplTime;
	SetStateSpaceModel();	// 状態空間モデルにセット
}

//! @brief 2慣性共振系のパラメータを設定する関数(Ds有り版)
//! @param[in]	Params		2慣性共振系パラメータ構造体
//! @param[in]	SmplTime	サンプリング周期 [s]
void TwoInertiaSimulator::SetParameters(const struct TwoInertiaParams& Params, const double SmplTime){
	Jl = Params.Jl;
	Dl = Params.Dl;
	Ds = Params.Ds;
	Ks = Params.Ks;
	Jm = Params.Jm;
	Dm = Params.Dm;
	Rg = Params.Rg;
	Kt = Params.Kt;
	Ts = SmplTime;
	SetStateSpaceModel();	// 状態空間モデルにセット
}

//! @brief モータ電流と負荷トルクを設定し状態を更新する関数
//! @param[in]	Current		q軸電流 [A]
//! @param[in]	LoadTorque	負荷トルク [Nm]
void TwoInertiaSimulator::SetCurrentAndLoadTorque(const double Current, const double LoadTorque){
	PlantSys.SetInput(Current, LoadTorque);	// 入力ベクトルの設定
	PlantSys.Update();						// 状態ベクトルを更新
	thl =  LoadInteg.GetSignal(PlantSys.GetOutput(1));	// 負荷側速度を積分して負荷側位置を計算
	thm = MotorInteg.GetSignal(PlantSys.GetOutput(3));	// モータ側速度を積分してモータ側位置を計算
}

//! @brief 現在の負荷側速度＆ねじれ角＆モータ側速度を取得する関数(引数で返す版)
//! @param[out]	LoadSpeed		負荷側速度 [rad/s]
//! @param[out]	TorsionAngle	ねじれ角 [rad]
//! @param[out]	MotorSpeed		モータ側速度 [rad/s]
void TwoInertiaSimulator::GetResponses(double& LoadSpeed, double& TorsionAngle, double& MotorSpeed){
	// 出力ベクトルから抽出して出力
	LoadSpeed    = PlantSys.GetOutput(1);
	TorsionAngle = PlantSys.GetOutput(2);
	MotorSpeed   = PlantSys.GetOutput(3);
}

//! @brief 現在の負荷側速度＆ねじれ角＆モータ側速度を取得する関数(タプルで返す版)
//! @return	タプル( 負荷側速度 [rad/s], ねじれ角 [rad], モータ側速度 [rad/s] )
std::tuple<double,double,double> TwoInertiaSimulator::GetResponses(void){
	// 出力ベクトルから抽出して出力
	const double LoadSpeed    = PlantSys.GetOutput(1);
	const double TorsionAngle = PlantSys.GetOutput(2);
	const double MotorSpeed   = PlantSys.GetOutput(3);
	return {LoadSpeed, TorsionAngle, MotorSpeed};
}

//! @brief 現在の負荷側位置を取得する関数
//! @return	負荷側位置応答 [rad]
double TwoInertiaSimulator::GetLoadPosition(void){
	return thl;
}

//! @brief 現在のモータ側位置を取得する関数
//! @return	モータ側位置応答 [rad]
double TwoInertiaSimulator::GetMotorPosition(void){
	return thm;
}

//! @brief 現在のモータ側速度を取得する関数
//! @return	モータ側速度応答 [rad/s]
double TwoInertiaSimulator::GetMotorSpeed(void){
	return PlantSys.GetOutput(3);
}

//! @brief 現在のねじれトルクを取得する関数
//! @return	ねじれトルク応答 [Nm]
double TwoInertiaSimulator::GetTorsionTorque(void){
	// 出力ベクトルから抽出
	const double LoadSpeed    = PlantSys.GetOutput(1);
	const double TorsionAngle = PlantSys.GetOutput(2);
	const double MotorSpeed   = PlantSys.GetOutput(3);
	
	// ねじれトルクを計算
	const double TorsionTorque = -Ds*LoadSpeed + Ks*TorsionAngle + Ds/Rg*MotorSpeed;
	return TorsionTorque;
}

//! @brief 電流と負荷を設定 → 状態更新 → 負荷側速度＆ねじれ角＆モータ側速度を取得する関数(引数で返す版)
//! @param[in]	Current			q軸電流 [A]
//! @param[in]	LoadTorque		負荷トルク [Nm]
//! @param[out]	LoadSpeed		負荷側速度 [rad/s]
//! @param[out]	TorsionAngle	ねじれ角 [rad]
//! @param[out]	MotorSpeed		モータ側速度 [rad/s]
void TwoInertiaSimulator::GetResponses(
	const double Current, const double LoadTorque, double& LoadSpeed, double& TorsionAngle, double& MotorSpeed){
	SetCurrentAndLoadTorque(Current, LoadTorque);		// 入力を設定して状態を更新
	GetResponses(LoadSpeed, TorsionAngle, MotorSpeed);	// 出力を計算
}

//! @brief 電流と負荷を設定 → 状態更新 → 負荷側速度＆ねじれ角＆モータ側速度を取得する関数(タプルで返す版)
//! @param[in]	Current		q軸電流 [A]
//! @param[in]	LoadTorque	負荷トルク [Nm]
//! @return	タプル( 負荷側速度 [rad/s], ねじれ角 [rad], モータ側速度 [rad/s] )
std::tuple<double,double,double> TwoInertiaSimulator::GetResponses(const double Current, const double LoadTorque){
	double LoadSpeed, TorsionAngle, MotorSpeed;
	GetResponses(Current, LoadTorque, LoadSpeed, TorsionAngle, MotorSpeed);
	return {LoadSpeed, TorsionAngle, MotorSpeed};
}

//! @brief 負荷側慣性を設定する関数
//! @param[in]	Inertia	負荷側慣性 [kgm^2]
void TwoInertiaSimulator::SetLoadInertia(const double Inertia){
	Jl = Inertia;
	SetStateSpaceModel();	// 状態空間モデルにセット
}

//! @brief ばね定数を設定する関数
//! @param[in]	Stiffness ばね定数 [Nm/rad]
void TwoInertiaSimulator::SetSpring(const double Stiffness){
	Ks = Stiffness;
	SetStateSpaceModel();	// 状態空間モデルにセット
}

//! @brief 状態空間モデルを設定する関数
void TwoInertiaSimulator::SetStateSpaceModel(void){
	// 連続系A行列の設定
	// 状態ベクトルの定義 xp = [ωl, θs, ωm]^T
	const Matrix<3,3> A = {
		-(Dl + Ds)/Jl,	Ks/Jl,			Ds/(Jl*Rg)			  ,
		-1,				0,				1.0/Rg				  ,
		Ds/(Jm*Rg),		-Ks/(Jm*Rg),	-(Ds/(Rg*Rg) + Dm)/Jm
	};
	
	// 連続系B行列の設定
	const Matrix<2,3> B = {
		0		, -1.0/Jl,
		0		, 0	     ,
		Kt/Jm	, 0
	};
	
	const auto C = Matrix<3,3>::eye();		// C行列の設定
	PlantSys.SetContinuous(A, B, C, Ts);	// 状態空間モデルの設定＆離散化
	LoadInteg.SetSmplTime(Ts);				// 積分器にサンプリング周期を設定
	MotorInteg.SetSmplTime(Ts);				// 積分器にサンプリング周期を設定
}

//! @brief シミュレータをリセットする関数
void TwoInertiaSimulator::Reset(void){
	PlantSys.ClearStateVector();
}

