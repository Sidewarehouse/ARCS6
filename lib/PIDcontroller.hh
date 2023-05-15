//! @file PIDcontroller.hh
//! @brief PID制御器クラス
//!
//! PID制御器 G(s) = Kp + Ki/s + Kd*s*gdis/(s+gdis) (双一次変換)
//!
//! @date 2019/07/29
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PIDCTRLLER
#define PIDCTRLLER

namespace ARCS {	// ARCS名前空間
	//! @brief PID制御器クラス
	class PIDcontroller {
		public:
			PIDcontroller(double Pgain, double Igain, double Dgain, double Bandwidth, double SmplTime);	//!< コンストラクタ
			PIDcontroller(PIDcontroller&& r);	//!< ムーブコンストラクタ
			~PIDcontroller();					//!< デストラクタ
			double GetSignal(double u);			//!< 出力信号の取得 u；入力信号
			void SetPgain(double Pgain);		//!< 比例ゲインの再設定 Pgain；比例ゲイン
			void SetIgain(double Igain);		//!< 積分ゲインの再設定 Igain；積分ゲイン
			void SetDgain(double Dgain);		//!< 微分ゲインの再設定 Dgain；微分ゲイン
			void SetBandwidth(double Bandwidth);//!< 擬似微分の帯域の再設定 Bandwidth；[rad/s] 帯域
			void SetSmplTime(double SmplTime);	//!< 制御周期の再設定 SmplTime；[s] 制御周期
			void ClearStateVars(void);			//!< すべての状態変数のリセット
			
		private:
			PIDcontroller(const PIDcontroller&) = delete;					//!< コピーコンストラクタ使用禁止
			const PIDcontroller& operator=(const PIDcontroller&) = delete;	//!< 代入演算子使用禁止
			double Ts;	//!< [s]		制御周期
			double Kp;	//!<			比例ゲイン
			double Ki;	//!<			積分ゲイン
			double Kd;	//!<			微分ゲイン
			double gpd;	//!< [rad/s]	擬似微分器の帯域
			double uZ1;	//!<			状態変数1 変数名Z1の意味はz変換のz^(-1)を示す
			double uZ2;	//!<			状態変数2
			double yZ1;	//!<			状態変数3
			double yZ2;	//!<			状態変数4
	};
}

#endif



