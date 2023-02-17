//! @file LoadVelocityObsrv.hh
//! @brief 負荷側速度オブザーバクラス
//!
//! モータ側速度とねじれトルクから，2慣性系の負荷側速度を推定する状態オブザーバ
//!
//! @date 2022/03/07
//! @author Kawai, Yusuke, and Yokokura, Yuki
//
// Copyright (C) 2011-2022 Kawai, Yusuke, and Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef LOADVELOCITYOBSRV
#define LOADVELOCITYOBSRV

#include <cassert>
#include <tuple>
#include "Matrix.hh"
#include "Discret.hh"
#include "StateSpaceSystem.hh"

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

namespace ARCS {	// ARCS名前空間
//! @brief 負荷側速度オブザーバクラス
template <size_t N = 1>
class LoadVelocityObsrv {
	public:
		//! @brief コンストラクタ
		//! @param[in]	GearRatio	減速比 [-]
		//! @param[in]	TorsionStiff	ねじれ剛性 [Nm/rad]
		//! @param[in]	Bandwidth	オブザーバの帯域 [rad/s]
		//! @param[in]	SmplTime	サンプリング周期 [s]
		LoadVelocityObsrv(const double GearRatio, const double TorsionStiff, const double Bandwidth, const double SmplTime)
			: Ks(TorsionStiff),
			  Rg(GearRatio),
			  g(Bandwidth),
			  Ts(SmplTime),
			  ObsrvSys()
		{
			// オブザーバゲイン
			const double k1 = Bandwidth;	// [rad/s]
			const double k2 = Bandwidth;	// [rad/s]
			const double k3 = Bandwidth;	// [rad/s]
			
			// 連続系オブザーバのA行列
			const Matrix<3,3> Ao = {
				-(k1 + k2 + k3)             , -Ks	 ,  0,
				(k1*k2 + k2*k3 + k3*k1)/Ks	, 0      ,  1,
				k1*k2*k3/Ks              	, 0      ,  0
			};
			
			// 連続系オブザーバのB行列
			const Matrix<2,3> Bo = {
				Ks/Rg,  ( k1 + k2 + k3 ),
				0     ,  -( k1*k2 + k2*k3 + k3*k1 )/Ks,
				0     ,  -k1*k2*k3/Ks
			};
			
			// オブザーバのC行列
			const Matrix<3,1> co = {
				0,  1,  0
			};
			
			// 状態空間モデルに設定
			ObsrvSys.SetContinuous(Ao, Bo, co, Ts);
			
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		LoadVelocityObsrv(LoadVelocityObsrv&& r)
			: Ks(r.Ks),
			  Rg(r.Rg),
			  g(r.g),
			  Ts(r.Ts),
			  ObsrvSys(r.ObsrvSys)
		{
			
		}
		
		//! @brief デストラクタ
		~LoadVelocityObsrv(){
			PassedLog();
		}
		
		//! @brief 推定負荷側速度を取得する関数
		//! @param[in]	MotorVelocity	モータ速度 [rad/s]
		//! @param[in]	TorsionTorque	ねじれトルク [Nm]
		//! @return	推定した負荷側速度 [rad/s]
		double GetLoadVelocity(const double MotorVelocity, const double TorsionTorque){
			const Matrix<1,2> u = {MotorVelocity, TorsionTorque};
			const Matrix<1,1> y = ObsrvSys.GetNextResponses(u);
			return y[1];
		}
		
	private:
		LoadVelocityObsrv(const LoadVelocityObsrv&) = delete;					//!< コピーコンストラクタ使用禁止
		const LoadVelocityObsrv& operator=(const LoadVelocityObsrv&) = delete;//!< 代入演算子使用禁止
		double Ks;	//!< [Nm/rad] 2慣性間の剛性
		double Rg;	//!< [-] 減速比
		double g;	//!< [rad/s] オブザーバの帯域
		double Ts;	//!< [s] サンプリング周期
		StateSpaceSystem<3,2,1> ObsrvSys;	//!< オブザーバの状態空間モデル
};
}

#endif

