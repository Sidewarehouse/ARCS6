//! @file TwoInertiaStateDistObsrv.hh
//! @brief 2慣性系状態外乱オブザーバクラス
//!
//! q軸電流とモータ側速度から，2慣性系の負荷側速度・ねじれ角・モータ側速度・負荷側ステップ外乱を推定する状態外乱オブザーバ
//! (普通の状態オブザーバと負荷側の0次外乱オブザーバを組み合わせたオブザーバ)
//!
//! @date 2020/04/06
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TWOINERTIASTATEDISTOBSRV
#define TWOINERTIASTATEDISTOBSRV

#include <cassert>
#include <tuple>
#include <ArcsMatrix.hh>
#include <Discret.hh>
#include "TwoInertiaParamDef.hh"
#include "Observer.hh"

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
//! @brief 2慣性系状態外乱オブザーバクラス
template <size_t N = 1>
class TwoInertiaStateDistObsrv {
	public:
		//! @brief コンストラクタ
		//! @param[in]	Params	2慣性共振系パラメータ構造体
		//! @param[in]	Bandwidth	オブザーバの帯域 [rad/s]
		//! @param[in]	SmplTime	サンプリング周期 [s]
		TwoInertiaStateDistObsrv(const struct TwoInertiaParamDef& Params, const double Bandwidth, const double SmplTime)
			: Kt(Params.Kt),
			  Jm(Params.Jm),
			  Dm(Params.Dm),
			  Jl(Params.Jl),
			  Dl(Params.Dl),
			  Ks(Params.Ks),
			  Rg(Params.Rg),
			  g(Bandwidth),
			  Ts(SmplTime),
			  A(),
			  b(),
			  c(),
			  k(),
			  StateDistObsrv()
		{
			CalcPlantStateModel();	// プラント状態空間モデルの計算
			CalcObserverGain();		// オブザーバゲインの計算
			StateDistObsrv.SetPlantModelAndGain(A, b, c, k, Ts);	// オブザーバの状態空間モデルの計算
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		TwoInertiaStateDistObsrv(TwoInertiaStateDistObsrv&& r)
			: Kt(r.Kt),
			  Jm(r.Jm),
			  Dm(r.Dm),
			  Jl(r.Jl),
			  Dl(r.Dl),
			  Ks(r.Ks),
			  Rg(r.Rg),
			  g(r.g),
			  Ts(r.Ts),
			  A(r.A),
			  b(r.b),
			  c(r.c),
			  k(r.k),
			  StateDistObsrv(std::move(r.StateDistObsrv))
		{
			
		}
		
		//! @brief デストラクタ
		~TwoInertiaStateDistObsrv(){
			PassedLog();
		}
		
		//! @brief 推定状態ベクトルを取得する関数(普通版)
		//! @param[in]	Current	q軸電流 [A]
		//! @param[in]	Velocity	モータ速度 [rad/s]
		//! @param[out]	xhat	状態ベクトル(負荷側外乱を含む) [ωl θs ωm τl]^T
		void GetEstimatedVect(const double Current, const double Velocity, ArcsMat<4,1>& xhat){
			const ArcsMat<2,1> u = {Current, Velocity};
			StateDistObsrv.Estimate(u, xhat);
		}
		
		//! @brief 推定状態ベクトルを取得する関数(ベクトルで返す版)
		//! @param[in]	Current	q軸電流 [A]
		//! @param[in]	Velocity	モータ速度 [rad/s]
		//! @return	推定した状態変数ベクトル(負荷側外乱を含む) [ωl θs ωm τl]^T
		ArcsMat<4,1> GetEstimatedVect(const double Current, const double Velocity){
			const ArcsMat<2,1> u = {Current, Velocity};
			return StateDistObsrv.Estimate(u);
		}
		
		//! @brief 各々の推定状態変数を取得する関数(タプル版)
		//! @param[in]	Current	q軸電流 [A]
		//! @param[in]	Velocity	モータ速度 [rad/s]
		//! @return	タプル( 負荷側速度 [rad/s], ねじれ角 [rad], モータ側速度 [rad/s], 負荷側外乱 [Nm] )
		std::tuple<double,double,double,double> GetEstimatedVars(const double Current, const double Velocity){
			const ArcsMat<2,1> u = {Current, Velocity};
			ArcsMat<4,1> xhat = StateDistObsrv.Estimate(u);
			return std::forward_as_tuple(xhat[1], xhat[2], xhat[3], xhat[4]);
		}
		
		//! @brief 推定状態ベクトルと負荷側外乱を分けて取得する関数(普通版)
		//! @param[in]	Current	q軸電流 [A]
		//! @param[in]	Velocity	モータ速度 [rad/s]
		//! @param[out]	xhat	状態ベクトル(負荷側外乱を含まない) [ωl θs ωm]^T
		//! @param[out]	LoadDisturbance	負荷側外乱トルク [Nm]
		void GetEstimatedVect(const double Current, const double Velocity, ArcsMat<3,1>& xhat, double& LoadDisturbance){
			const ArcsMat<2,1> u = {Current, Velocity};
			ArcsMat<4,1> xhat_dis;
			StateDistObsrv.Estimate(u, xhat_dis);
			xhat[1] = xhat_dis[1];
			xhat[2] = xhat_dis[2];
			xhat[3] = xhat_dis[3];
			LoadDisturbance = xhat_dis[4];
		}
		
	private:
		TwoInertiaStateDistObsrv(const TwoInertiaStateDistObsrv&) = delete;					//!< コピーコンストラクタ使用禁止
		const TwoInertiaStateDistObsrv& operator=(const TwoInertiaStateDistObsrv&) = delete;//!< 代入演算子使用禁止
		double Kt;	//!< [Nm/A] トルク定数
		double Jm;	//!< [kgm^2] モータ慣性
		double Dm;	//!< [Nm/(rad/s)] モータ粘性
		double Jl;	//!< [kgm^2] 負荷側慣性
		double Dl;	//!< [Nm/(rad/s)] 負荷側粘性
		double Ks;	//!< [Nm/rad] 2慣性間の剛性
		double Rg;	//!< [-] 減速比
		double g;	//!< [rad/s] オブザーバの帯域
		double Ts;	//!< [s] サンプリング周期
		ArcsMat<4,4> A;	//!< プラントのA行列
		ArcsMat<4,1> b;	//!< プラントのBベクトル
		ArcsMat<1,4> c;	//!< プラントのCベクトル
		ArcsMat<4,1> k;	//!< オブザーバのゲインベクトル
		Observer<4> StateDistObsrv;	//!< 状態外乱オブザーバ
		
		//! @brief プラントの状態空間モデルを計算する関数
		void CalcPlantStateModel(void){
			// プラントの状態方程式 (x = [ωl θs ωm τl]^T, u = iqref, y = ωm)
			A.Set(
				-Dl/Jl,  Ks/Jl		,  0		, -1.0/Jl,
				-1    ,	 0			,  1.0/Rg	,  0     ,
				 0	  ,	-Ks/(Jm*Rg)	, -Dm/Jm	,  0     ,
				 0	  ,  0			,  0		,  0	 
			);
			b.Set(
				0,
				0,
				Kt/Jm,
				0
			);
			c.Set(0, 0, 1, 0);
		}
		
		//! @brief オブザーバゲインを計算する関数
		void CalcObserverGain(void){
			// 推定帯域gで重根配置されるオブザーバゲインを計算
			const double k1 = -(Jm*Rg*(Dl - 2.0*Jl*g)*(Dl*Dl - 2.0*Dl*Jl*g + 2.0*Jl*Jl*g*g - 2.0*Ks*Jl))/(Jl*Jl*Jl*Ks);
			const double k2 = (- Jm*Dl*Dl*Rg*Rg + 4.0*Jm*Dl*Jl*Rg*Rg*g - 6.0*Jm*Jl*Jl*Rg*Rg*g*g + Ks*Jl*Jl + Jm*Ks*Jl*Rg*Rg)/(Jl*Jl*Ks*Rg);
			const double k3 = -(Dl*Jm + Dm*Jl - 4.0*Jl*Jm*g)/(Jl*Jm);
			const double k4 = -(Jl*Jm*Rg*g*g*g*g)/Ks;
			
			// オブザーバゲインベクトルにセット
			k.Set(
				k1,
				k2,
				k3,
				k4
			);
		}
};
}

#endif

