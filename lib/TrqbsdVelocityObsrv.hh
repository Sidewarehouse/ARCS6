//! @file TrqbsdVelocityObsrv.hh
//! @brief トルクセンサベース速度オブザーバ
//! @date 2018/12/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2018 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TRQBSDVELOCITYOBSRV
#define TRQBSDVELOCITYOBSRV

#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief トルクセンサベース速度オブザーバ
	class TrqbsdVelocityObsrv {
		public:
			TrqbsdVelocityObsrv(
				double TrqConst, double Inertia, double Viscosity, double GearRatio, double Bandwidth, double SmplTime
			);	//!< コンストラクタ
			~TrqbsdVelocityObsrv();	//!< デストラクタ
			double GetVelocity(double Current, double TorsionTorque, double Position);	//!< 速度推定値を計算する関数
			
		private:
			TrqbsdVelocityObsrv(const TrqbsdVelocityObsrv&) = delete;					//!< コピーコンストラクタ使用禁止
			const TrqbsdVelocityObsrv& operator=(const TrqbsdVelocityObsrv&) = delete;	//!< 代入演算子使用禁止
			
			const double Ktn;	//!< [Nm/A]	トルク定数
			const double Jmn;	//!< [kgm^2]	モータ側慣性モーメント
			const double Dmn;	//!< [Nm/(rad/s)]	モータ側粘性
			const double Rg;	//!< [-]		減速比
			const double w;		//!< [rad/s]	推定帯域
			const double Ts;	//!< [s]		サンプリング時間
			const double e;		//!< e = exp(-Dmn/Jmn*Ts)
			const double k1;	//!< オブザーバゲイン1
			const double k2;	//!< オブザーバゲイン2
			Matrix<1,2> k;		//!< オブザーバゲインベクトル
			Matrix<2,2> Ad;		//!< 離散系 A行列
			Matrix<2,2> Bd;		//!< 離散系 B行列
			Matrix<2,1> c;		//!< C行列
			Matrix<1,2> u;		//!< プラント入力ベクトル
			Matrix<1,2> xh;		//!< 状態推定ベクトル
			Matrix<1,2> xh_next;//!< 次の時刻の状態推定ベクトル
			Matrix<1,1> y;		//!< プラント出力ベクトル
	};
}

#endif

