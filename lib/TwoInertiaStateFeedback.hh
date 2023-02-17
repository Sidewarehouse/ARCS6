//! @file TwoInertiaStateFeedback.hh
//! @brief 2慣性共振系の状態フィードバックレギュレータ
//!
//! 2慣性共振系に対して，状態オブザーバの推定値を用いて状態フィードバックを掛けて任意の極配置にセットするレギュレータ。
//! (MATLABでいうところの「reg」のようなもの)
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TWOINERTIASTATEFEEDBACK
#define TWOINERTIASTATEFEEDBACK

#include <cassert>
#include <complex>
#include "Matrix.hh"
#include "TwoInertiaParamDef.hh"
#include "TwoInertiaStateObsrv.hh"

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
//! @brief 2慣性共振系の状態フィードバックレギュレータ
//! @tparam 
//template <>
class TwoInertiaStateFeedback {
	public:
		//! @brief コンストラクタ(実数極版)
		//! @param[in]	Params	2慣性共振系パラメータ構造体
		//! @param[in]	Pole1,2,3	状態フィードバックの指定極(1～3 実数極) [rad/s]
		//! @param[in]	Bandwidth	状態オブザーバの帯域 [rad/s]
		//! @param[in]	SmplTime	サンプリング周期 [s]
		TwoInertiaStateFeedback(
			const struct TwoInertiaParamDef& Params,
			const double Pole1, const double Pole2, const double Pole3,
			const double Bandwidth, const double SmplTime
		)
			: Kt(Params.Kt),
			  Jm(Params.Jm),
			  Dm(Params.Dm),
			  Jl(Params.Jl),
			  Dl(Params.Dl),
			  Ks(Params.Ks),
			  Rg(Params.Rg),
			  gsob(Bandwidth),
			  Ts(SmplTime),
			  f(),
			  SOB(Params, Bandwidth, SmplTime)
		{
			CalcStateFeedbackGain(Pole1, Pole2, Pole3);
			PassedLog();
		}
		
		//! @brief コンストラクタ(実数極版)
		//! @param[in]	Params	2慣性共振系パラメータ構造体
		//! @param[in]	Pole1	状態フィードバックの指定極1(実数極) [rad/s]
		//! @param[in]	Pole2	状態フィードバックの指定極2(複素極) [rad/s]
		//! @param[in]	Pole3	状態フィードバックの指定極3(複素極) [rad/s]
		//! @param[in]	Bandwidth	状態オブザーバの帯域 [rad/s]
		//! @param[in]	SmplTime	サンプリング周期 [s]
		TwoInertiaStateFeedback(
			const struct TwoInertiaParamDef& Params,
			const double Pole1, const std::complex<double>& Pole2, const std::complex<double>& Pole3,
			const double Bandwidth, const double SmplTime
		)
			: Kt(Params.Kt),
			  Jm(Params.Jm),
			  Dm(Params.Dm),
			  Jl(Params.Jl),
			  Dl(Params.Dl),
			  Ks(Params.Ks),
			  Rg(Params.Rg),
			  gsob(Bandwidth),
			  Ts(SmplTime),
			  f(),
			  SOB(Params, Bandwidth, SmplTime)
		{
			CalcStateFeedbackGain(Pole1, Pole2, Pole3);
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		TwoInertiaStateFeedback(TwoInertiaStateFeedback&& r)
			: Kt(r.Kt),
			  Jm(r.Jm),
			  Dm(r.Dm),
			  Jl(r.Jl),
			  Dl(r.Dl),
			  Ks(r.Ks),
			  Rg(r.Rg),
			  gsob(r.gsob),
			  Ts(r.Ts),
			  f(r.f),
			  SOB(std::move(r.SOB))
		{
			
		}

		//! @brief デストラクタ
		~TwoInertiaStateFeedback(){
			PassedLog();
		}
		
		//! @brief 状態フィードバックレギュレータの補償電流を取得する関数
		//! @param[in]	Current	q軸電流指令 [A]
		//! @param[in]	Velocity	モータ側速度 [rad/s]
		double GetFeedbackCurrent(const double Current, const double Velocity){
			const Matrix<1,1> ireg = f*SOB.GetEstimatedVect(Current, Velocity);
			return ireg[1];
		}
		
	private:
		TwoInertiaStateFeedback(const TwoInertiaStateFeedback&) = delete;					//!< コピーコンストラクタ使用禁止
		const TwoInertiaStateFeedback& operator=(const TwoInertiaStateFeedback&) = delete;	//!< 代入演算子使用禁止
		double Kt;	//!< [Nm/A] トルク定数
		double Jm;	//!< [kgm^2] モータ慣性
		double Dm;	//!< [Nm/(rad/s)] モータ粘性
		double Jl;	//!< [kgm^2] 負荷側慣性
		double Dl;	//!< [Nm/(rad/s)] 負荷側粘性
		double Ks;	//!< [Nm/rad] 2慣性間の剛性
		double Rg;	//!< [-] 減速比
		double gsob;//!< [rad/s] 状態オブザーバの帯域
		double Ts;	//!< [s] サンプリング周期
		Matrix<3,1> f;				//!< 状態フィードバックゲインベクトル
		TwoInertiaStateObsrv<> SOB;	//!< 状態オブザーバ
		
		//! @brief 状態フィードバックゲインを計算する関数(実数極版)
		//! @param[in]	p1	指定極1(実数) [rad/s]
		//! @param[in]	p2	指定極2(実数) [rad/s]
		//! @param[in]	p3	指定極3(実数) [rad/s]
		void CalcStateFeedbackGain(double p1, double p2, double p3){
			const double fl = -(Jm*Rg*(Dl*Dl*Dl - 2.0*Dl*Jl*Ks + Dl*Dl*Jl*p1 + Dl*Dl*Jl*p2 + Dl*Dl*Jl*p3 - Jl*Jl*Ks*p1 - Jl*Jl*Ks*p2 - Jl*Jl*Ks*p3 + Dl*Jl*Jl*p1*p2 + Dl*Jl*Jl*p1*p3 + Dl*Jl*Jl*p2*p3 + Jl*Jl*Jl*p1*p2*p3))/(Jl*Jl*Ks*Kt);
			const double fs =  (Jm*Rg*(p1*p2 + p1*p3 + p2*p3 - Ks/Jl + Dl*Dl/(Jl*Jl) - Ks/(Jm*Rg*Rg) + (Dl*p1)/Jl + (Dl*p2)/Jl + (Dl*p3)/Jl))/Kt;
			const double fm = -(Dl*Jm + Dm*Jl + Jl*Jm*p1 + Jl*Jm*p2 + Jl*Jm*p3)/(Jl*Kt);
			f.Set(fl, fs, fm);
		}
		
		//! @brief 状態フィードバックゲインを計算する関数(複素極版)
		//! @param[in]	p1	指定極1(実数) [rad/s]
		//! @param[in]	p2	指定極2(複素数) [rad/s]
		//! @param[in]	p3	指定極3(複素数) [rad/s]
		void CalcStateFeedbackGain(double p1, const std::complex<double>& p2, const std::complex<double>& p3){
			arcs_assert(p2 == conj(p3));	// 極がペアになってるかチェック
			const std::complex<double> fl = -(Jm*Rg*(Dl*Dl*Dl - 2.0*Dl*Jl*Ks + Dl*Dl*Jl*p1 + Dl*Dl*Jl*p2 + Dl*Dl*Jl*p3 - Jl*Jl*Ks*p1 - Jl*Jl*Ks*p2 - Jl*Jl*Ks*p3 + Dl*Jl*Jl*p1*p2 + Dl*Jl*Jl*p1*p3 + Dl*Jl*Jl*p2*p3 + Jl*Jl*Jl*p1*p2*p3))/(Jl*Jl*Ks*Kt);
			const std::complex<double> fs =  (Jm*Rg*(p1*p2 + p1*p3 + p2*p3 - Ks/Jl + Dl*Dl/(Jl*Jl) - Ks/(Jm*Rg*Rg) + (Dl*p1)/Jl + (Dl*p2)/Jl + (Dl*p3)/Jl))/Kt;
			const std::complex<double> fm = -(Dl*Jm + Dm*Jl + Jl*Jm*p1 + Jl*Jm*p2 + Jl*Jm*p3)/(Jl*Kt);
			f.Set(fl.real(), fs.real(), fm.real());
		}
		
};
}

#endif

