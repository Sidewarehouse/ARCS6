//! @file Integrator2.hh
//! @brief 2次積分器クラス
//!
//! 2次積分器 G(s)=1/s^2 (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef INTEGRATOR2
#define INTEGRATOR2

namespace ARCS {	// ARCS名前空間
	//! @brief 2次積分器クラス
	class Integrator2 {
		private:
			double Ts;	// [s]	制御周期
			double uZ1;	// 		状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double uZ2;	// 		状態変数2 変数名Z2の意味はz変換のz^(-2)を示す
			double yZ1;	//		状態変数3
			double yZ2;	//		状態変数4
		
		public:
			explicit Integrator2(double SmplTime);	// コンストラクタ SmplTime；[s] 制御周期
			~Integrator2();						// デストラクタ
			double GetSignal(double u);			// 出力信号の取得 u；入力信号
			void SetSmplTime(double SmplTime);	// 制御周期の再設定 SmplTime；[s] 制御周期
			void ClearStateVars(void);			// すべての状態変数のリセット
	};
}

#endif



