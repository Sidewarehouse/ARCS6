//! @file P_Dcontroller.hh
//! @brief P-D制御器クラス(微分先行型)
//!
//! 微分先行型の比例-微分制御器
//!
//! @date 2022/03/30
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef P_DCONTROLLER
#define P_DCONTROLLER

#include <cassert>
#include "MovingDifferentiator.hh"

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
//! @brief P-D制御器クラス(微分先行型)
//! @tparam	W	微分の窓幅(差分取るときのサンプル数) [-] (デフォルト値 = 4)
template <size_t W = 4>
class P_Dcontroller {
	public:
		//! @brief 空コンストラクタ
		P_Dcontroller(void)
			: Diff(), Kp(0), Kd(0)
		{
			
		}
		
		//! @brief コンストラクタ
		//! @param[in]	Pgain	比例ゲイン [*]
		//! @param[in]	Dgain	微分ゲイン [*]
		P_Dcontroller(const double Pgain, const double Dgain)
			: Diff(), Kp(Pgain), Kd(Dgain)
		{
			EventLogVar(Kp);
			EventLogVar(Kd);
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		P_Dcontroller(P_Dcontroller&& r)
			 : Diff(r.Diff), Kp(r.Kp), Kd(r.Kd)
		{
			
		}

		//! @brief デストラクタ
		~P_Dcontroller(){
			
		}
		
		//! @brief P-D制御器の出力信号を取得する関数
		//! y = (ref - res)*Kp - d(res)/d(time)*Kd
		//! @param[in]	ref		指令値 [*]
		//! @param[in]	res		応答値 [*]
		//! @param[in]	time	現在時刻 [s] (微分計算を正確にするために必要)
		//! @return	P-D制御器の出力値 [*]
		double GetSignal(const double ref, const double res, const double time){
			return (ref - res)*Kp - Diff.GetSignal(res, time)*Kd;
		}
		
		//! @brief Pゲインを設定する関数
		//! @param[in]	Pgain	比例ゲイン [*]
		void SetPgain(const double Pgain){
			Kp = Pgain;
		}
		
		//! @brief Dゲインを設定する関数
		//! @param[in]	Dgain	微分ゲイン [*]
		void SetDgain(const double Dgain){
			Kd = Dgain;
		}
		
		//! @brief PゲインとDゲインを設定する関数
		//! @param[in]	Pgain	比例ゲイン [*]
		//! @param[in]	Dgain	微分ゲイン [*]
		void SetPDgain(const double Pgain, const double Dgain){
			SetPgain(Pgain);
			SetDgain(Dgain);
		}
		
	private:
		P_Dcontroller(const P_Dcontroller&) = delete;					//!< コピーコンストラクタ使用禁止
		const P_Dcontroller& operator=(const P_Dcontroller&) = delete;	//!< 代入演算子使用禁止
		MovingDifferentiator<W> Diff;	//!< 微分器
		double Kp;	//!< 比例ゲイン
		double Kd;	//!< 微分ゲイン
};
}

#endif

