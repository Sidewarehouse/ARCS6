//! @file TwoInertiaSimulator.hh
//! @brief 2慣性共振系シミュレータV2(ベクトル版)
//!
//! 複数の2慣性共振系を模擬して，負荷側速度，ねじれ角，モータ速度を計算して出力する。
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TWOINERTIASIMULATORS
#define TWOINERTIASIMULATORS

#include <cassert>
#include <array>
#include "Matrix.hh"
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

namespace ARCS {	// ARCS名前空間
//! @brief 2慣性共振系シミュレータV2(ベクトル版)
//! @tparam N	軸数
template <size_t N>
class TwoInertiaSimulators {
	public:
		//! @brief コンストラクタ
		//! @param[in]	Params	2慣性系パラメータ構造体の配列
		//! @param[in]	SmplTime	サンプリング周期 [s]
		TwoInertiaSimulators(const std::array<struct TwoInertiaParamDef, N>& Params, const double SmplTime)
			: PlantSyss()
		{
			PassedLog();
			// 各軸数分の2慣性系パラメータ設定
			for(size_t i = 0; i < N; ++i){
				PlantSyss.at(i).SetParameters(Params.at(i), SmplTime);
			}
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		TwoInertiaSimulators(TwoInertiaSimulators&& r)
			: PlantSyss(std::move(r.PlantSyss))
		{
			
		}

		//! @brief デストラクタ
		~TwoInertiaSimulators(){
			PassedLog();
		}
		
		//! @brief 予めモータ電流と負荷トルクを設定する関数(std::array版)
		//! @param[in]	Current	q軸電流の配列 [A]
		//! @param[in]	LoadTorque	負荷トルクの配列 [Nm]
		void SetCurrentAndLoadTorque(const std::array<double, N>& Current, const std::array<double, N>& LoadTorque){
			// 各軸の電流と負荷トルクを設定
			for(size_t i = 0; i < N; ++i){
				PlantSyss.at(i).SetCurrentAndLoadTorque(Current.at(i), LoadTorque.at(i));
			}
		}
		
		//! @brief 予めモータ電流と負荷トルクを設定する関数(ベクトル版)
		//! @param[in]	Current	q軸電流の縦ベクトル [A]
		//! @param[in]	LoadTorque	負荷トルクの縦ベクトル [Nm]
		void SetCurrentAndLoadTorque(const Matrix<1,N>& Current, const Matrix<1,N>& LoadTorque){
			std::array<double, N> iq, ts;
			Current.StoreArray(iq);
			LoadTorque.StoreArray(ts);
			SetCurrentAndLoadTorque(iq, ts);
		}
		
		//! @brief 予め設定された電流と負荷から負荷側速度＆ねじれ角＆モータ側速度を計算する関数(std::array版)
		//! @param[out]	LoadSpeed	負荷側速度の配列 [rad/s]
		//! @param[out]	TorsionAngle	ねじれ角の配列 [rad]
		//! @param[out]	MotorSpeed	モータ側速度の配列 [rad/s]
		void GetResponses(std::array<double, N>& LoadSpeed, std::array<double, N>& TorsionAngle, std::array<double, N>& MotorSpeed){
			// 各軸の応答を計算
			for(size_t i = 0; i < N; ++i){
				PlantSyss.at(i).GetResponses(LoadSpeed.at(i), TorsionAngle.at(i), MotorSpeed.at(i));
			}
		}
		
		//! @brief 予め設定された電流と負荷から負荷側速度＆ねじれ角＆モータ側速度を計算する関数(ベクトル版)
		//! @param[out]	LoadSpeed	負荷側速度の縦ベクトル [rad/s]
		//! @param[out]	TorsionAngle	ねじれ角の縦ベクトル [rad]
		//! @param[out]	MotorSpeed	モータ側速度の縦ベクトル [rad/s]
		void GetResponses(Matrix<1,N>& LoadSpeed, Matrix<1,N>& TorsionAngle, Matrix<1,N>& MotorSpeed){
			std::array<double, N> wl, ths, wm;
			GetResponses(wl, ths, wm);
			LoadSpeed.LoadArray(wl);
			TorsionAngle.LoadArray(ths);
			MotorSpeed.LoadArray(wm);
		}
		
		//! @brief モータ側位置を取得する関数(std::array版)
		//! 注意： GetResponses()を事前に実行しないと位置は反映されない。
		//! @param[out]	MotorPosition	モータ側位置応答の配列 [rad]
		void GetMotorPosition(std::array<double, N>& MotorPosition){
			// 各軸の位置を取得
			for(size_t i = 0; i < N; ++i){
				MotorPosition.at(i) = PlantSyss.at(i).GetMotorPosition();
			}
		}
		
		//! @brief モータ側位置を取得する関数(ベクトル版)
		//! 注意： GetResponses()を事前に実行しないと位置は反映されない。
		//! @param[out]	MotorPosition	モータ側位置応答の縦ベクトル [rad]
		void GetMotorPosition(Matrix<1,N>& MotorPosition){
			std::array<double, N> thm;
			GetMotorPosition(thm);
			MotorPosition.LoadArray(thm);
		}
		
		//! @brief 負荷側位置を取得する関数(std::array版)
		//! 注意： GetResponses()を事前に実行しないと位置は反映されない。
		//! @param[out]	LoadPosition	負荷側位置応答の配列 [rad]
		void GetLoadPosition(std::array<double, N>& LoadPosition){
			// 各軸の位置を取得
			for(size_t i = 0; i < N; ++i){
				LoadPosition.at(i) = PlantSyss.at(i).GetLoadPosition();
			}
		}
		
		//! @brief 負荷側位置を取得する関数(ベクトル版)
		//! 注意： GetResponses()を事前に実行しないと位置は反映されない。
		//! @param[out]	LoadPosition	負荷側位置応答の縦ベクトル [rad]
		void GetLoadPosition(Matrix<1,N>& LoadPosition){
			std::array<double, N> thl;
			GetLoadPosition(thl);
			LoadPosition.LoadArray(thl);
		}
		
	private:
		TwoInertiaSimulators(const TwoInertiaSimulators&) = delete;					//!< コピーコンストラクタ使用禁止
		const TwoInertiaSimulators& operator=(const TwoInertiaSimulators&) = delete;//!< 代入演算子使用禁止
		
		std::array<TwoInertiaSimulator, N> PlantSyss;	//!< 2慣性系シミュレータの配列
		
};
}

#endif

