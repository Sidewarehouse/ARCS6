//! @file LoadsideDistObsrv.hh
//! @brief 同一次元負荷側外乱オブザーバ
//!
//! モータ側速度と捻れトルクから負荷側外乱トルクを推定します。
//!
//! @date 2018/12/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2018 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef LOADSIDEDISTOBSRV
#define LOADSIDEDISTOBSRV

#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief 同一次元負荷側外乱オブザーバ
	class LoadsideDistObsrv {
		public:
			LoadsideDistObsrv(double Stiffness, double LoadInertia, double GearRatio, double Bandwidth, double SmplTime);	//!< コンストラクタ(オブザーバの極を3重根にする場合)
			~LoadsideDistObsrv();	//!< デストラクタ
			double GetLoadsideTorque(double MotorSpeed, double TorsionTorque);	//!< 負荷側トルクを推定する関数
			
		private:
			LoadsideDistObsrv(const LoadsideDistObsrv&) = delete;					// コピーコンストラクタ使用禁止
			const LoadsideDistObsrv& operator=(const LoadsideDistObsrv&) = delete;	// 代入演算子使用禁止
			double Ksn;	//!< [Nm/rad] 2慣性間の剛性
			double Jln;	//!< [kgm^2] 負荷側慣性
			double Rg;	//!< [-] 減速比
			double l1;	//!< [rad/s] オブザーバの極1
			double l2;	//!< [rad/s] オブザーバの極2
			double l3;	//!< [rad/s] オブザーバの極3
			double Ts;	//!< [s] 制御周期
			
			Matrix<3,3> Ac;	//!<  連続系 A行列
			Matrix<2,3> Bc;	//!<  連続系 B行列
			Matrix<3,3> Ad;	//!<  離散系 A行列
			Matrix<2,3> Bd;	//!<  離散系 B行列
			Matrix<3,1> c;	//!<  C行列
			Matrix<1,2> u;	//!<  入力ベクトル
			Matrix<1,3> x;	//!<  状態ベクトル
			Matrix<1,3> x_next;	//!<  次の時刻の状態ベクトル
			Matrix<1,1> y;	//!<  出力ベクトル
	};
}

#endif

