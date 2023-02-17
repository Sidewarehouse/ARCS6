//! @file ResoLoopCutDistObsrv.hh
//! @brief 共振ループ切断外乱オブザーバ(ベクトル版)
//!
//! モータ側トルク指令とねじれトルクセンサ値から，モータ側外乱と負荷側速度外乱を推定するオブザーバで，
//! 2慣性共振系の共振ループを切断するための特殊な外乱オブザーバのクラス
//!
//! @date 2022/03/30
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef RESOLOOPCUTDISTOBSRVS
#define RESOLOOPCUTDISTOBSRVS

#include <cassert>
#include <array>
#include "ResoLoopCutDistObsrv.hh"

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
//! @brief 共振ループ切断外乱オブザーバ(ベクトル版)
//! @tparam	N	軸数
template <size_t N>
class ResoLoopCutDistObsrvs {
	public:
		//! @brief コンストラクタ
		//! @param[in]	Params		2慣性共振系パラメータ構造体の配列
		//! @param[in]	Bandwidth	推定帯域 [rad/s]
		//! @param[in]	SmplTime	サンプリング周期 [s]
		ResoLoopCutDistObsrvs(const std::array<struct TwoInertiaParams, N>& Params, const double Bandwidth, const double SmplTime)
			: RLCDObs()
		{
			for(size_t i = 0; i < N; ++i) RLCDObs.at(i).SetParameters(Params.at(i), Bandwidth, SmplTime);
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		ResoLoopCutDistObsrvs(ResoLoopCutDistObsrvs&& r)
			: RLCDObs(r.RLCDObs)
		{
			
		}
		
		//! @brief デストラクタ
		~ResoLoopCutDistObsrvs(){
			PassedLog();
		}
		
		//! @brief 推定帯域を設定する関数
		//! @param[in]	Bandwidth	推定帯域 [rad/s]
		void SetBandwidth(const double Bandwidth){
			for(size_t i = 0; i < N; ++i) RLCDObs.at(i).SetBandwidth(Bandwidth);
		}
		
		//! @brief モータ側トルク補償値を取得する関数(引数で返す版)
		//! @param[in]	TorqueRef		モータ側トルク指令 [Nm]
		//! @param[in]	TorsionTorque	ねじれトルク [Nm]
		//! @param[out]	モータ側トルク補償値 [Nm]
		void GetCompTorque(const Matrix<1,N>& TorqueRef, const Matrix<1,N>& TorsionTorque, Matrix<1,N>& CompTorque){
			for(size_t i = 1; i <= N; ++i) CompTorque[i] = RLCDObs.at(i - 1).GetCompTorque(TorqueRef[i], TorsionTorque[i]);
		}
		
		//! @brief モータ側トルク補償値を取得する関数(戻り値で返す版)
		//! @param[in]	TorqueRef		モータ側トルク指令 [Nm]
		//! @param[in]	TorsionTorque	ねじれトルク [Nm]
		//! @param[out]	モータ側トルク補償値 [Nm]
		Matrix<1,N> GetCompTorque(const Matrix<1,N>& TorqueRef, const Matrix<1,N>& TorsionTorque){
			Matrix<1,N> ret;
			GetCompTorque(TorqueRef, TorsionTorque, ret);
			return ret;
		}
		
		//! @brief RLC-TTC用のP-D制御ゲインを計算する関数(仮設)
		//! @param[in]	wt	ねじれトルク制御帯域 [rad/s]
		//! @param[in]	zt	ねじれトルク制御制動係数 [-]
		//! @param[out]	Kpt	P-Dねじれトルク制御器 比例ゲイン
		//! @param[out] Kdt P-Dねじれトルク制御器 微分ゲイン
		void GetPDgainForRLCTTC(const double wt, const double zt, Matrix<1,N>& Kpt, Matrix<1,N>& Kdt){
			for(size_t i = 1; i <= N; ++i) RLCDObs.at(i - 1).GetPDgainForRLCTTC(wt, zt, Kpt[i], Kdt[i]);
		}
		
	private:
		ResoLoopCutDistObsrvs(const ResoLoopCutDistObsrvs&) = delete;					//!< コピーコンストラクタ使用禁止
		const ResoLoopCutDistObsrvs& operator=(const ResoLoopCutDistObsrvs&) = delete;	//!< 代入演算子使用禁止
		std::array<ResoLoopCutDistObsrv, N> RLCDObs;									//!< 共振ループ切断外乱オブザーバの配列
};
}

#endif

