//! @file I-P-I-Pcontroller.cc
//! @brief I-P-I-P制御器
//! @date 2019/09/19
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef IPIPCONTROLLER
#define IPIPCONTROLLER

#include "Integrator.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief I-P-I-P制御器
	//! @tparam	T	積分方法の指定
	template <IntegralType T>
	class IPIPcontroller {
		public:
			//! @brief コンストラクタ
			//! @param[in]	Igain1	積分ゲイン1
			//! @param[in]	Pgain1	比例ゲイン1
			//! @param[in]	Igain2	積分ゲイン2
			//! @param[in]	Pgain2	比例ゲイン2
			//! @param[in]	SmplTime	制御周期 [s]
			IPIPcontroller(double Igain1, double Pgain1, double Igain2, double Pgain2, double SmplTime)
				: Ki1(Igain1), Kp1(Pgain1), Ki2(Igain2), Kp2(Pgain2), Ts(SmplTime), Integ1(SmplTime), Integ2(SmplTime)
			{
				
			}
			
			//! @brief ムーブコンストラクタ
			//! @param[in]	r	右辺値
			IPIPcontroller(IPIPcontroller&& r)
				: Ki1(r.Ki1), Kp1(r.Kp1), Ki2(r.Ki2), Kp2(r.Kp2), Ts(r.Ts), Integ1(r.Integ1), Integ2(r.Integ2)
			{
				
			}
			
			//! @brief デストラクタ
			~IPIPcontroller(){
				
			}

			//! @brief I-P-I-P制御器出力信号の計算
			//! @param[in]	ref	指令値
			//! @param[in]	res	応答値
			//! @return	制御器出力
			double GetSignal(double ref, double res){
				return Integ2.GetSignal(  Integ1.GetSignal(ref - res)*Ki1 - Kp1*res  )*Ki2 - Kp2*res;
			}
		
		private:
			IPIPcontroller(const IPIPcontroller&) = delete;						//!< コピーコンストラクタ使用禁止
			const IPIPcontroller& operator=(const IPIPcontroller&) = delete;	//!< 代入演算子使用禁止
			double Ki1;				//!< 積分ゲイン1
			double Kp1;				//!< 比例ゲイン1
			double Ki2;				//!< 積分ゲイン2
			double Kp2;				//!< 比例ゲイン2
			double Ts;				//!< 制御周期
			Integrator<T> Integ1;	//!< 積分器1
			Integrator<T> Integ2;	//!< 積分器2
	};
}

#endif

