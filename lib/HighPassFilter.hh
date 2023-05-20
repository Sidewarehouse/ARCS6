//! @file HighPassFilter.hh
//! @brief 高域通過濾波器クラス
//!
//! 1次高域通過濾波器 G(s)=s/(s+g) (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef HIGHPASSFILTER
#define HIGHPASSFILTER

namespace ARCS {	// ARCS名前空間
	//! @brief 高域通過濾波器クラス
	class HighPassFilter {
		private:
			double Ts;	// [s]		制御周期
			double g;	// [rad/s]	遮断周波数
			double uZ1;	// 			状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double yZ1;	//			状態変数2
		
		public:
			HighPassFilter(double CutFreq, double SmplTime);
			// コンストラクタ CutFreq；[rad/s] 遮断周波数，SmplTime；[s] 制御周期
			~HighPassFilter();					// デストラクタ
			double GetSignal(double u);			// 出力信号の取得 u；入力信号
			void SetCutFreq(double CutFreq);	// 遮断周波数の再設定 CutFreq；[rad/s] 遮断周波数
			void SetSmplTime(double SmplTime);	// 制御周期の再設定 SmplTime；[s] 制御周期
			void ClearStateVars(void);			// すべての状態変数のリセット
	};
}

#endif



