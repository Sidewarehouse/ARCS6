//! @file HighPassFilter2.hh
//! @brief 2次高域通過濾波器クラス
//!
//! 2次高域通過濾波器 G(s)=s^2/(s^2 + w/Q*s + w^2) (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef HIGHPASSFILTER2
#define HIGHPASSFILTER2

namespace ARCS {	// ARCS名前空間
	//! @brief 2次高域通過濾波器クラス
	class HighPassFilter2 {
		private:
			double Ts;	// [s]		制御周期
			double w;	// [rad/s]	遮断周波数
			double Q;	//			鋭さ(0.5のときに重根設定となる)
			double uZ1;	//			状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double uZ2;	//			状態変数2
			double uZ3;	//			状態変数3
			double yZ1;	//			状態変数4
			double yZ2;	//			状態変数5
			double yZ3;	//			状態変数6
		
		public:
			HighPassFilter2(double CutFreq, double Sharp, double SmplTime);
			// コンストラクタ CutFreq；[rad/s] 遮断周波数，Sharp；鋭さ，SmplTime；[s] 制御周期
			~HighPassFilter2();					// デストラクタ
			double GetSignal(double u);			// 出力信号の取得 u；入力信号
			void SetCutFreq(double CutFreq);	// 遮断周波数の再設定 CutFreq；[rad/s] 遮断周波数
			void SetSharpness(double Sharp);	// 鋭さの再設定 Sharp；鋭さ
			void SetSmplTime(double SmplTime);	// 制御周期の再設定 SmplTime；[s] 制御周期
			void ClearStateVars(void);			// すべての状態変数のリセット
	};
}

#endif



