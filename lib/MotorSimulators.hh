//! @file MotorSimulator.hh
//! @brief モータシミュレータ(ベクトル対応版)
//!
//! 複数のモータを模擬する
//!
//! @date 2020/04/14
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef MOTORSIMULATORS
#define MOTORSIMULATORS

#include <cassert>
#include <array>
#include "Matrix.hh"
#include "MotorParamDef.hh"
#include "MotorSimulator.hh"

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
//! @brief モータシミュレータ(ベクトル対応版)
//! @tparam N	軸数
template <size_t N>
class MotorSimulators {
	public:
		//! @brief コンストラクタ
		//! @param[in]	Params	モータパラメータ構造体の配列
		//! @param[in]	SmplTime	サンプリング周期 [s]
		MotorSimulators(const std::array<struct MotorParamDef, N>& Params, const double SmplTime)
			: PlantSyss()
		{
			PassedLog();
			// 各軸数分のモータシミュレータのパラメータ設定
			for(size_t i = 0; i < N; ++i){
				PlantSyss.at(i).SetParameters(Params.at(i).Kt, Params.at(i).Jm, Params.at(i).Dm, SmplTime);
			}
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		MotorSimulators(MotorSimulators&& r)
			: PlantSyss(std::move(r.PlantSyss))
		{
			
		}

		//! @brief デストラクタ
		~MotorSimulators(){
			PassedLog();
		}
		
		//! @brief 予め電流と負荷トルクを設定する関数
		//! @param[in]	Current	[A] 電流の配列
		//! @param[in]	LoadTorque	[Nm] 負荷トルクの配列
		void SetCurrentAndLoadTorque(const std::array<double, N> Current, const std::array<double, N> LoadTorque){
			// 各軸の電流と負荷トルクを設定
			for(size_t i = 0; i < N; ++i){
				PlantSyss.at(i).SetCurrentAndLoadTorque(Current.at(i), LoadTorque.at(i));
			}
		}
		
		//! @brief 予めセットされたモータ電流と負荷トルクから速度と位置を計算する関数
		//! @param[out]	Velocity	[rad/s] 速度
		//! @param[out] Position	[rad] 位置
		void GetVelocityAndPosition(std::array<double, N>& Velocity, std::array<double, N>& Position){
			// 各軸の応答計算
			for(size_t i = 0; i < N; ++i){
				PlantSyss.at(i).GetVelocityAndPosition(Velocity.at(i), Position.at(i));
			}
		}
		
		
	private:
		MotorSimulators(const MotorSimulators&) = delete;					//!< コピーコンストラクタ使用禁止
		const MotorSimulators& operator=(const MotorSimulators&) = delete;	//!< 代入演算子使用禁止
		
		std::array<MotorSimulator, N> PlantSyss;	//!< モータシミュレータの配列
};
}

#endif

