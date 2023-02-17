//! @file PDcontroller.hh
//! @brief PD制御器クラス
//!
//! PD制御器 G(s) = Kp + Kd*s*gdis/(s+gdis) (双一次変換)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PDCTRLLER
#define PDCTRLLER

namespace ARCS {	// ARCS名前空間
	//! @brief PD制御器クラス
	class PDcontroller {
		private:
			double Ts;	// [s]		制御周期
			double Kp;	//			比例ゲイン
			double Kd;	//			微分ゲイン
			double gpd;	// [rad/s]	擬似微分器の帯域
			double uZ1;	// 			状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double yZ1;	//			状態変数2
		
		public:
			PDcontroller(double Pgain, double Dgain, double Bandwidth, double SmplTime);
			// コンストラクタ
			// Pgain；比例ゲイン，Dgain；微分ゲイン，Bandwidth；[rad/s] 帯域，SmplTime；[s] 制御周期
			~PDcontroller();					// デストラクタ
			double GetSignal(double u);			// 出力信号の取得 u；入力信号
			void SetPgain(double Pgain);		// 比例ゲインの再設定 Pgain；比例ゲイン
			void SetDgain(double Dgain);		// 微分ゲインの再設定 Dgain；微分ゲイン
			void SetBandwidth(double Bandwidth);// 擬似微分の帯域の再設定 Bandwidth；[rad/s] 帯域
			void SetSmplTime(double SmplTime);	// 制御周期の再設定 SmplTime；[s] 制御周期
			void ClearStateVars(void);			// すべての状態変数のリセット
	};
}

#endif



