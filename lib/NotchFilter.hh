//! @file NotchFilter.hh
//! @brief ノッチフィルタクラス
//!
//! ノッチフィルタ G(s)=( s^2 + w^2 )/( s^2 + w/Q*s + w^2) (双一次変換)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef NOTCHFILTER
#define NOTCHFILTER

namespace ARCS {	// ARCS名前空間
	//! @brief ノッチフィルタクラス
	class NotchFilter {
		private:
			double Ts;	// [s]		制御周期
			double w;	// [rad/s]	遮断中心周波数
			double Q;	//			鋭さ
			double uZ1;	//			状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double uZ2;	//			状態変数2
			double yZ1;	//			状態変数3
			double yZ2;	//			状態変数4
		
		public:
			NotchFilter(double AngFreq, double Sharp, double SmplTime);
			// コンストラクタ  AngFreq；[rad/s] 遮断中心周波数，Sharp；鋭さ，SmplTime；[s] 制御周期
			~NotchFilter();						// デストラクタ
			double GetSignal(double u);			// 出力信号の取得 u；入力信号
			void SetCuttoff(double AngFreq);	// 遮断中心周波数の再設定 AngFreq；[rad/s] 帯域
			void SetSharpness(double Sharp);	// 鋭さの再設定 Sharp；鋭さ
			void SetSmplTime(double SmplTime);	// 制御周期の再設定 SmplTime；[s] 制御周期
			void ClearStateVars(void);			// すべての状態変数のリセット
	};
}

#endif



