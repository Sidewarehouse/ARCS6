//! @file P_Dcontrollers.hh
//! @brief P-D制御器クラス(微分先行型)(ベクトル版)
//!
//! 微分先行型の比例-微分制御器 多軸対応バージョン
//!
//! @date 2022/03/30
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef P_DCONTROLLERS
#define P_DCONTROLLERS

#include <cassert>
#include <array>
#include "Matrix.hh"
#include "P-Dcontroller.hh"

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
//! @tparam	N	軸数
//! @tparam W	微分の窓幅(差分取るときのサンプル数) [-] (デフォルト値 = 4)
template <size_t N, size_t W = 4>
class P_Dcontrollers {
	public:
		//! @brief 空コンストラクタ
		P_Dcontrollers(void)
			: P_Dcon()
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ(std::array版)
		//! @param[in]	Pgain	比例ゲイン [*]
		//! @param[in]	Dgain	微分ゲイン [*]
		P_Dcontrollers(const std::array<double, N>& Pgain, const std::array<double, N>& Dgain)
			: P_Dcon()
		{
			SetPgain(Pgain);
			SetDgain(Dgain);
			PassedLog();
		}

		//! @brief コンストラクタ(Matrix版)
		//! @param[in]	Pgain	比例ゲイン [*]
		//! @param[in]	Dgain	微分ゲイン [*]
		P_Dcontrollers(const Matrix<1,N>& Pgain, const Matrix<1,N>& Dgain)
			: P_Dcon()
		{
			SetPgain(Pgain);
			SetDgain(Dgain);
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		P_Dcontrollers(P_Dcontrollers&& r)
			 : P_Dcon(r.P_Dcon)
		{
			
		}

		//! @brief デストラクタ
		~P_Dcontrollers(){
			PassedLog();
		}
		
		//! @brief P-D制御器の出力信号を取得する関数(引数で返す版)
		//! out = (ref - res)*Kp - d(res)/d(time)*Kd
		//! @param[in]	ref		指令値 [*]
		//! @param[in]	res		応答値 [*]
		//! @param[in]	time	現在時刻 [s] (微分計算を正確にするために必要)
		//! @param[out]	out		P-D制御器の出力値 [*]
		void GetSignal(const Matrix<1,N>& ref, const Matrix<1,N>& res, const double time, Matrix<1,N>& out){
			for(size_t i = 1; i <= N; ++i){
				out[i] = P_Dcon.at(i-1).GetSignal(ref[i], res[i], time);
			}
		}
		
		//! @brief P-D制御器の出力信号を取得する関数(返り値で返す版)
		//! out = (ref - res)*Kp - d(res)/d(time)*Kd
		//! @param[in]	ref		指令値 [*]
		//! @param[in]	res		応答値 [*]
		//! @param[in]	time	現在時刻 [s] (微分計算を正確にするために必要)
		//! @return	P-D制御器の出力値(縦ベクトル) [*]
		Matrix<1,N> GetSignal(const Matrix<1,N>& ref, const Matrix<1,N>& res, const double time){
			Matrix<1,N> ret;
			GetSignal(ref, res, time, ret);
			return ret;
		}
		
		//! @brief Pゲインを設定する関数(std::array版)
		//! @param[in]	Pgain	比例ゲイン [*]
		void SetPgain(const std::array<double, N>& Pgain){
			for(size_t i = 0; i < N; ++i) P_Dcon.at(i).SetPgain(Pgain.at(i));	// 各々のPゲインに設定
		}
		
		//! @brief Dゲインを設定する関数(std::array版)
		//! @param[in]	Dgain	微分ゲイン [*]
		void SetDgain(const std::array<double, N>& Dgain){
			for(size_t i = 0; i < N; ++i) P_Dcon.at(i).SetDgain(Dgain.at(i));	// 各々のDゲインに設定
		}
		
		//! @brief Pゲインを設定する関数(Matrix版)
		//! @param[in]	Pgain	比例ゲイン [*]
		void SetPgain(const Matrix<1,N>& Pgain){
			for(size_t i = 1; i <= N; ++i) P_Dcon.at(i - 1).SetPgain(Pgain[i]);	// 各々のPゲインに設定
		}
		
		//! @brief Dゲインを設定する関数(Matrix版)
		//! @param[in]	Dgain	微分ゲイン [*]
		void SetDgain(const Matrix<1,N>& Dgain){
			for(size_t i = 1; i <= N; ++i) P_Dcon.at(i - 1).SetDgain(Dgain[i]);	// 各々のDゲインに設定
		}
		
		//! @brief PDゲインを設定する関数(std::array版)
		//! @param[in]	Pgain	比例ゲイン [*]
		//! @param[in]	Dgain	微分ゲイン [*]
		void SetPDgain(const std::array<double, N>& Pgain, const std::array<double, N>& Dgain){
			for(size_t i = 0; i < N; ++i) P_Dcon.at(i).SetPDgain(Pgain.at(i), Dgain.at(i));	// 各々のPDゲインに設定
		}
		
		//! @brief PDゲインを設定する関数(Matrix版)
		//! @param[in]	Pgain	比例ゲイン [*]
		//! @param[in]	Dgain	微分ゲイン [*]
		void SetPDgain(const Matrix<1,N>& Pgain, const Matrix<1,N>& Dgain){
			for(size_t i = 1; i <= N; ++i) P_Dcon.at(i - 1).SetPDgain(Pgain[i], Dgain[i]);	// 各々のPDゲインに設定
		}
		
	private:
		P_Dcontrollers(const P_Dcontrollers&) = delete;					//!< コピーコンストラクタ使用禁止
		const P_Dcontrollers& operator=(const P_Dcontrollers&) = delete;//!< 代入演算子使用禁止
		std::array<P_Dcontroller<W>, N> P_Dcon;							//!< P-D制御器の配列
};
}

#endif

