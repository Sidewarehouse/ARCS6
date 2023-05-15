//! @file MovingAverage.hh
//! @brief 移動平均クラス
//!
//! 入力信号の移動平均を計算し出力するクラス。
//!
//! @date 2020/04/15
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef MOVINGAVERAGE
#define MOVINGAVERAGE

#include <cassert>
#include "RingBuffer.hh"
#include "Statistics.hh"

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
//! @brief 移動平均クラス
//! @tparam N	移動平均する時刻方向のデータ数
template <size_t N>
class MovingAverage {
	public:
		//! @brief コンストラクタ
		MovingAverage()
			: WindowData()
		{
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		MovingAverage(MovingAverage&& r)
			: WindowData(std::move(r.WindowData))
		{
			
		}

		//! @brief デストラクタ
		~MovingAverage(){
			PassedLog();
		}
		
		//! @brief 移動平均の先頭データに入力信号をセットして出力信号を得る関数
		//! @param[in]	u	入力信号
		//! @return	出力信号
		double GetSignal(const double u){
			WindowData.SetFirstValue(u);				// リングバッファの先頭に入力信号を詰める
			return Statistics::Mean(WindowData.Buffer);	// リングバッファ内の平均値を返す
		}
		
		//! @brief 移動平均の窓データを指定した値で埋める関数
		//! @param[in]	u	指定値
		void Fill(const double u){
			WindowData.FillBuffer(u);	// リングバッファを指定した値で埋める
		}
		
		//! @brief 移動平均の窓データをクリアする関数
		void Clear(void){
			WindowData.ClearBuffer();
		}
		
	private:
		MovingAverage(const MovingAverage&) = delete;					//!< コピーコンストラクタ使用禁止
		const MovingAverage& operator=(const MovingAverage&) = delete;	//!< 代入演算子使用禁止
		
		RingBuffer<double, N, false> WindowData;	//!< 移動する窓データ
};
}

#endif

