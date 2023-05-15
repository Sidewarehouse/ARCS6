//! @file PIcontroller.cc
//! @brief PI制御器(通常版，積分停止法版，リミット偏差フィードバック版)
//! @date 2019/09/20
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PICONTROLLER
#define PICONTROLLER

#include "Integrator.hh"
#include "Limiter.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief PI制御器
	//! @tparam	T	積分方法の指定(デフォルトは後退オイラー)
	template <IntegralType T = IntegralType::BACKWARD_EULER>
	class PIcontroller {
		public:
			//! @brief コンストラクタ
			//! @param[in]	Pgain	比例ゲイン
			//! @param[in]	Igain	積分ゲイン
			//! @param[in]	SmplTime	制御周期 [s]
			PIcontroller(const double Pgain, const double Igain, const double SmplTime)
				: Kp(Pgain), Ki(Igain), Integ(SmplTime)
			{
				
			}
			
			//! @brief ムーブコンストラクタ
			//! @param[in]	r	右辺値
			PIcontroller(PIcontroller&& r)
				: Kp(r.Kp), Ki(r.Ki), Integ(r.Integ)
			{
				
			}
			
			//! @brief デストラクタ
			~PIcontroller(){
				
			}

			//! @brief PI制御器出力信号の計算(通常版)
			//! @param[in]	u	制御器入力
			//! @return	制御器出力
			double GetSignal(const double u){
				return Kp*u + Ki*Integ.GetSignal(u);
			}
			
			//! @brief PI制御器出力信号の計算(積分停止法版)
			//! @param[in]	u	制御器入力
			//! @param[in]	Limit	積分停止リミット値
			//! @return	制御器出力
			double GetSignal(const double u, const double Limit){
				Integ.Enable(false);						// 積分器を一旦停止
				double ytest = Kp*u + Ki*Integ.GetSignal(u);// 制御器出力を積分値を更新せずに一旦計算
				if(-Limit < ytest && ytest < Limit){		// リミットの範囲内だったら，
					Integ.Enable(true);						// 積分器を再開して，
					return Kp*u + Ki*Integ.GetSignal(u);	// 普通にPIの計算して返す
				}
				return Limiter(ytest, Limit);				// リミットの外だったらリミット値を返す
			}
			
			//! @brief PI制御器出力信号の計算(リミット偏差フィードバック版)
			//! @param[in]	u		制御器入力
			//! @param[in]	LimitIn		リミット入力側
			//! @param[in]	LimitOut	リミット出力側
			//! @return	制御器出力
			double GetSignal(const double u, const double LimitIn, const double LimitOut){
				return Kp*u + Ki*Integ.GetSignal( u - 1.0/Kp*(LimitIn - LimitOut) );
			}
			
			//! @brief 積分器を動かすか止めるかの選択
			//! @param[in] Flag trueのとき積分器作動，falseのとき積分停止
			void Enable(const bool Flag){
				Integ.Enable(Flag);
			}
			
			//! @brief 積分器のリセット
			void ClearIntegralValue(void){
				Integ.ClearIntegralValue();
			}
			
		private:
			PIcontroller(const PIcontroller&) = delete;						//!< コピーコンストラクタ使用禁止
			const PIcontroller& operator=(const PIcontroller&) = delete;	//!< 代入演算子使用禁止
			double Kp;				//!< 比例ゲイン
			double Ki;				//!< 積分ゲイン
			Integrator<T> Integ;	//!< 積分器
	};
}

#endif

