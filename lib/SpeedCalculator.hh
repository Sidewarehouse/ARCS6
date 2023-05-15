//! @file SpeedCalculator.hh
//! @brief 速度計算器
//!
//! 位置計測値と時刻計測値から速度を計算する
//!
//! @date 2020/04/07
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SPEEDCALCULATOR
#define SPEEDCALCULATOR

#include <cmath>
#include "RingBuffer.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief 速度計算器
	//! @tparam	N	差分のサンプリング数 [-]
	template <size_t N>
	class SpeedCalculator {
		public:
			//! @brief コンストラクタ
			SpeedCalculator()
				: PositionBuff(), TimeBuff(), FirstTime(true)
			{
				
			}

			//! @brief ムーブコンストラクタ
			//! @param[in]	r	右辺値
			SpeedCalculator(SpeedCalculator&& r)
				: PositionBuff(r.PositionBuff), TimeBuff(r.TimeBuff), FirstTime(r.FirstTime)
			{
				
			}

			//! @brief デストラクタ
			~SpeedCalculator(){
				
			}
			
			//! @brief 速度を計算する関数
			//! @param[in]	Position	現在位置 [rad],[m]
			//! @param[in]	Time		現在時刻 [s]
			//! @return	速度 [rad/s], [m/s]
			double GetSpeed(const double Position, const double Time){
				if(FirstTime){
					// 初めて呼ばれた場合は，
					PositionBuff.FillBuffer(Position);	// 初期位置でバッファを埋めておく
					TimeBuff.FillBuffer(Time);			// 初期時刻でバッファを埋めておく
					FirstTime = false;					// フラグリセット
					return 0;							// 最初の1回目は速度を計算できないのでゼロ
				}else{
					// 2回目以降は，
					const double dx = Position - PositionBuff.GetFinalValue();	// 現在位置と位置リングバッファの最後尾との偏差の計算
					const double dt = Time - TimeBuff.GetFinalValue();			// 現在時刻と時刻リングバッファの最後尾との偏差の計算
					const double v = dx/dt;					// 速度の計算
					PositionBuff.SetFirstValue(Position);	// 位置リングバッファの先頭に現在位置を詰める
					TimeBuff.SetFirstValue(Time);			// 時刻リングバッファの先頭に現在時刻を詰める
					if(std::isfinite(v)){
						return v;	// 速度が有限値なら計算結果を返す
					}else{
						return 0;	// ゼロ割のときはゼロにしておく
					}
				}
			}
			
			//! @brief リセット
			void Reset(void){
				FirstTime = true;
			}
			
		private:
			SpeedCalculator(const SpeedCalculator&) = delete;					//!< コピーコンストラクタ使用禁止
			const SpeedCalculator& operator=(const SpeedCalculator&) = delete;	//!< 代入演算子使用禁止
			RingBuffer<double, N, false> PositionBuff;	//!< 位置データ用リングバッファ
			RingBuffer<double, N, false> TimeBuff;		//!< 時間データ用リングバッファ
			bool FirstTime;								//!< 初めて呼ばれたかどうかのフラグ
	};
}

#endif

