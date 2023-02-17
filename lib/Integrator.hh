//! @file Integrator.hh
//! @brief 積分器(テンプレート版)
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef INTEGRATOR
#define INTEGRATOR

#include <cassert>

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
//! @brief 積分方法の定義
enum class IntegralType {
	FOWARD_EULER,	//!< 前進オイラー(前進差分)
	BACKWARD_EULER,	//!< 後退オイラー(後退差分)
	TRAPEZOIDAL		//!< 台形積分(双1次変換, Tustin変換)
};

//! @brief 積分器
//! @tparam	T	積分方法の指定(デフォルトは後退オイラー)
template <IntegralType T = IntegralType::BACKWARD_EULER>
class Integrator {
	public:
		//! @brief 空コンストラクタ
		Integrator(void)
			: Ts(0),			// [s] 制御周期の格納
			  uZ1(0), yZ1(0),	// 状態変数の初期化
			  Ena(true)			// 積分器作動フラグ
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ
		//! @param[in]	SmplTime	制御周期 [s]
		explicit Integrator(const double SmplTime)
			: Ts(SmplTime),		// [s] 制御周期の格納
			  uZ1(0), yZ1(0),	// 状態変数の初期化
			  Ena(true)			// 積分器作動フラグ
		{
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		Integrator(Integrator&& r)
			: Ts(r.Ts),					// [s] 制御周期の格納
			  uZ1(r.uZ1), yZ1(r.yZ1),	// 状態変数の初期化
			  Ena(r.Ena)				// 積分器作動フラグ
		{
			
		}
		
		//! @brief デストラクタ
		~Integrator(){
			PassedLog();
		}

		//! @brief 出力信号の取得
		//! @param[in]	u	入力信号
		//! @return 積分器出力
		double GetSignal(const double u){
			double y = 0;
			
			if(Ena == false)return yZ1;	// 積分器停止指令の場合は積分動作を停止して積分値をそのまま出力
			
			// 積分方法によって計算を変える
			switch(T){
				case IntegralType::FOWARD_EULER:
					y = Ts*uZ1 + yZ1;			// 前進オイラー(前進差分)
					break;
				case IntegralType::BACKWARD_EULER:
					y = Ts*u + yZ1;				// 後退オイラー(後退差分)
					break;
				case IntegralType::TRAPEZOIDAL:
					y = Ts/2.0*(u + uZ1) + yZ1;	// 台形積分(双1次変換, Tustin変換)
					break;
				default:
					y = Ts*u + yZ1;				// デフォルトは後退オイラー
					break;
			}
			
			// 状態変数の更新
			uZ1 = u;
			yZ1 = y;
			
			return y;
		}

		//! @brief 制御周期の再設定
		//! @param[in] SmplTime 制御周期[s]
		void SetSmplTime(const double SmplTime){
			Ts = SmplTime;	// [s] 制御周期の再設定
		}

		//! @brief 積分初期値の設定
		//! @param[in] InitValue 初期値
		void SetInitial(const double InitValue){
			yZ1 = InitValue;
		}
		
		//! @brief 積分値の設定
		//! @param[in] Value 初期値
		void SetIntegralValue(const double Value){
			yZ1 = Value;
		}
		
		//! @brief 積分器を動かすか止めるかの選択
		//! @param[in] Flag trueのとき積分器作動，falseのとき積分停止
		void Enable(const bool Flag){
			Ena = Flag;
		}

		//! @brief 積分器のリセット
		void ClearIntegralValue(void){
			uZ1 = 0;	// 状態変数1のゼロクリア
			yZ1 = 0;	// 状態変数2のゼロクリア
		}
		
	private:
		Integrator(const Integrator&) = delete;					//!< コピーコンストラクタ使用禁止
		const Integrator& operator=(const Integrator&) = delete;//!< 代入演算子使用禁止
		double Ts;	//!< [s]制御周期
		double uZ1;	//!< 状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
		double yZ1;	//!< 状態変数2
		bool Ena;	//!< 積分器作動指令 trueのとき積分器作動，falseのとき積分停止
	
};
}

#endif

