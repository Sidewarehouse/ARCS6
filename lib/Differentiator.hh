//! @file Differentiator.hh
//! @brief  擬似微分器クラス
//!
//! 擬似微分器 G(s)=(s*gpd)/(s+gpd) (双一次変換)
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef DIFFERENTIATOR
#define DIFFERENTIATOR

namespace ARCS {	// ARCS名前空間
	//! @brief  擬似微分器クラス
	class Differentiator {
		private:
			double Ts;	// [s]		制御周期
			double gpd;	// [rad/s]	擬似微分の帯域
			double uZ1;	//			状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double yZ1;	//			状態変数2
		
		public:
			Differentiator(double Bandwidth, double SmplTime);
			// コンストラクタ Bandwidth；[rad/s] 帯域，SmplTime；[s] 制御周期
			~Differentiator();					// デストラクタ
			double GetSignal(double u);			// 出力信号の取得 u；入力信号
			void SetBandwidth(double Bandwidth);// 擬似微分の帯域の再設定 Bandwidth；[rad/s] 帯域
			void SetSmplTime(double SmplTime);	// 制御周期の再設定 SmplTime；[s] 制御周期
			void ClearStateVars(void);			// すべての状態変数のリセット
	};
}

#endif



