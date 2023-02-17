//! @file HysteresisComparator.hh
//! @brief ヒステリシス比較器
//! @date 2019/02/22
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef HYSTERESISCOMPARATOR
#define HYSTERESISCOMPARATOR

namespace ARCS {	// ARCS名前空間
	//! @brief ヒステリシス比較器
	class HysteresisComparator {
		public:
			HysteresisComparator(double Left, double Right);//!< コンストラクタ
			~HysteresisComparator();	//!< デストラクタ
			double GetSignal(double u);	//!< 出力信号の取得
			
		private:
			HysteresisComparator(const HysteresisComparator&) = delete;					// コピーコンストラクタ使用禁止
			const HysteresisComparator& operator=(const HysteresisComparator&) = delete;// 代入演算子使用禁止
			double R;	//!< ヒステリシス特性の右側の比較値
			double L;	//!< ヒステリシス特性の左側の比較値
			double yZ1;	//!< 状態変数
	};
}

#endif

