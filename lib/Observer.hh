//! @file Observer.hh
//! @brief 一般形のオブザーバクラス(とりあえずSISO版)
//!
//! 任意の制御対象の状態空間モデルからオブザーバを構成して，入出力信号から状態ベクトルを推定する。
//!
//! @date 2022/06/29
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef OBSERVER
#define OBSERVER

#include <cassert>
#include <complex>
#include "Matrix.hh"
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
//! @brief 一般形のオブザーバクラス
//! @tparam N	制御対象の次数
template <size_t N>
class Observer {
	public:
		//! @brief 空コンストラクタ
		Observer(void)
			: ObsrvSys()
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ
		Observer(const Matrix<N,N>& A, const Matrix<1,N>& b, const Matrix<N,1>& c)
			: ObsrvSys()
		{
			Matrix<1,N,std::complex<double>> lambda = eigen(A);
			PrintMat(tp(lambda));
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		Observer(Observer&& r)
				: ObsrvSys(std::move(r.ObsrvSys))
		{
			
		}

		//! @brief デストラクタ
		~Observer(){
			PassedLog();
		}
		
		//! @brief 対象のプラントモデルとゲインを設定する関数
		//! @param[in]	A	プラントの連続系A行列
		//! @param[in]	b	プラントの連続系bベクトル
		//! @param[in]	c	プラントのcベクトル
		//! @param[in]	k	オブザーバゲインベクトル
		//! @param[in]	Ts	サンプリング周期 [s]
		void SetPlantModelAndGain(const Matrix<N,N>& A, const Matrix<1,N>& b, const Matrix<N,1>& c, const Matrix<1,N>& k, const double Ts){
			const Matrix<N,N> Ao = A - k*c;			// オブザーバの連続系A行列
			Matrix<2,N> Bo;							// オブザーバの連続系B行列
			setcolumn(Bo, b, 1);					// オブザーバの連続系B行列1列目
			setcolumn(Bo, k, 2);					// オブザーバの連続系B行列2列目
			const auto Co = Matrix<N,N>::eye();		// オブザーバの出力行列
			ObsrvSys.SetContinuous(Ao, Bo, Co, Ts);	// オブザーバの状態空間モデルを設定
		}
		
		//! @brief 状態推定の計算をして状態ベクトルを返す関数(普通版)
		//! @param[in]	u	オブザーバの入力ベクトル
		//! @param[out]	xhat	推定状態ベクトル
		void Estimate(const Matrix<1,2>& u, Matrix<1,N>& xhat){
			ObsrvSys.GetResponses(u, xhat);			// 推定演算
		}
		
		//! @brief 状態推定の計算をして状態ベクトルを返す関数(ベクトルを返す版)
		//! @param[in]	u	オブザーバの入力ベクトル
		//! @return 推定状態ベクトル
		Matrix<1,N> Estimate(const Matrix<1,2>& u){
			return ObsrvSys.GetResponses(u);		// 推定演算
		}
		
	private:
		Observer(const Observer&) = delete;					//!< コピーコンストラクタ使用禁止
		const Observer& operator=(const Observer&) = delete;//!< 代入演算子使用禁止
		StateSpaceSystem<N,2,N> ObsrvSys;	//!< オブザーバの状態空間モデル
};
}

#endif

