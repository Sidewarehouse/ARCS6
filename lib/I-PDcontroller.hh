//! @file I-PDcontroller.hh
//! @brief I-PD制御器
//! @date 2020/03/05
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef I_PDCONTROLLER
#define I_PDCONTROLLER

#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief I-PD制御器
	class I_PDcontroller {
		public:
			I_PDcontroller(double Igain, double Pgain, double Dgain, double Bandwidth, double SmplTime);	//!< コンストラクタ
			~I_PDcontroller();	//!< デストラクタ
			double GetSignal(double ref, double res);//!< 制御器の出力信号を計算して出力する関数 ref : 指令値, res : 応答値
			void SetIgain(double Igain);	//!< 積分ゲインを再設定する関数
			void SetPgain(double Pgain);	//!< 先行比例ゲインを再設定する関数
			void SetDgain(double Dgain);	//!< 先行微分ゲインを再設定する関数
			void SetBandwidth(double Bandwidth);	//!< 擬似微分の帯域を再設定する関数
			
		private:
			I_PDcontroller(const I_PDcontroller&) = delete;					//!< コピーコンストラクタ使用禁止
			const I_PDcontroller& operator=(const I_PDcontroller&) = delete;//!< 代入演算子使用禁止
			double Ki;	//!< 積分ゲイン
			double Kp;	//!< 先行比例ゲイン
			double Kd;	//!< 先行微分ゲイン
			double g;	//!< [rad/s] 擬似微分器の帯域
			double Ts;	//!< [s] サンプリング周期
			Matrix<2,2> Ac;		//!< 連続系A行列
			Matrix<2,2> Bc;		//!< 連続系B行列
			Matrix<2,2> Ad;		//!< 離散系A行列
			Matrix<2,2> Bd;		//!< 離散系B行列
			Matrix<2,1> c;		//!< C行列
			Matrix<2,1> d;		//!< D行列
			Matrix<1,2> u;		//!< 入力ベクトル
			Matrix<1,2> x;		//!< 状態ベクトル
			Matrix<1,2> x_next;	//!< 次の時刻の状態ベクトル
			Matrix<1,1> y;		//!< 出力ベクトル
	};
}

#endif



