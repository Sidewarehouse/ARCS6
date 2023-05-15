//! @file FixedAverage.hh
//! @brief FixedAverageクラス
//!
//! 時刻ゼロから現在時刻までの信号の平均を計算する。
//!
//! @date 2019/07/29
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef FIXEDAVERAGE
#define FIXEDAVERAGE

namespace ARCS {	// ARCS名前空間
	//! @brief FixedAverageクラス
	class FixedAverage {
		public:
			FixedAverage();					//!< コンストラクタ
			FixedAverage(FixedAverage&& r);	//!< ムーブコンストラクタ
			~FixedAverage();				//!< デストラクタ
			double GetSignal(double in);	//!< ひたすら平均値を計算し続ける関数
			void Reset(void);				//!< 平均値とカウント数をリセットする関数
			
		private:
			FixedAverage(const FixedAverage&) = delete;					//!< コピーコンストラクタ使用禁止
			const FixedAverage& operator=(const FixedAverage&) = delete;//!< 代入演算子使用禁止
			unsigned int count;	//!< 関数が呼ばれた回数
			double ave;			//!< 時刻ゼロからの平均値
	};
}

#endif



