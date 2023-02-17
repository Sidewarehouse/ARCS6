//! @file PhaseLeadLag.hh
//! @brief 位相進み/遅れ補償器クラス
//!
//! 位相補償器 G(s)=(1+s*a/g)/(1+s/g) (双一次変換)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PHASELEADLAG
#define PHASELEADLAG

namespace ARCS {	// ARCS名前空間
	//! @brief 位相進み/遅れ補償器クラス
	class PhaseLeadLag {
		private:
			double Ts;	// [s]		制御周期
			double g;	// [rad/s]	ゲイン特性の折れ点周波数
			double a;	//			a < 1 で位相遅れ，1 < a で位相進み
			double uZ1;	// 			状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double yZ1;	//			状態変数2
		
		public:
			PhaseLeadLag(double AngFreq, double Phase, double SmplTime);
			// コンストラクタ
			// AngFreq；[rad/s] ゲイン特性の折れ点周波数，a < 1 で位相遅れ，1 < a で位相進み，SmplTime；[s] 制御周期
			~PhaseLeadLag();					// デストラクタ
			double GetSignal(double u);			// 出力信号の取得 u；入力信号
			void SetSmplTime(double SmplTime);	// 制御周期の再設定 SmplTime；[s] 制御周期
			void SetFreq(double AngFreq);			// ゲイン特性の折れ点周波数の再設定 AngFreq；[rad/s] 折れ点周波数
			void SetPhase(double Phase);		// a < 1 で位相遅れ，1 < a で位相進み
			void ClearStateVars(void);			// すべての状態変数のリセット
	};
}

#endif



