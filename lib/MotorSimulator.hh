//! @file MotorSimulator.hh
//! @brief モータシミュレータ
//!
//! モータを模擬する
//!
//! @date 2020/04/14
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef MOTORSIMULATOR
#define MOTORSIMULATOR

#include "Matrix.hh"
#include "StateSpaceSystem.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief モータシミュレータ
	class MotorSimulator {
		public:
			MotorSimulator(void);					//!< 空コンストラクタ
			MotorSimulator(const double TrqConst, const double MotorInert, const double MotorFric, const double SmplTime);	//!< コンストラクタ
			MotorSimulator(MotorSimulator&& r);		//!< ムーブコンストラクタ
			~MotorSimulator();						//!< デストラクタ
			void SetParameters(const double TrqConst, const double MotorInert, const double MotorFric, const double SmplTime); //!< パラメータを設定する関数
			void SetCurrentAndLoadTorque(const double current, const double loadtorque);							//!< 予めモータ電流と負荷トルクを設定する関数
			void GetResponses(const double current, const double loadtorque, double& velocity, double& position);	//!< モータ速度と位置を返す関数(普通版)
			std::tuple<double,double> GetResponses(const double current, const double loadtorque);	//!< モータ速度と位置を計算する関数(タプル版)
			void GetVelocityAndPosition(double& velocity, double& position);	//!< 予めセットされたモータ電流と負荷トルクから速度と位置を返す関数(普通版)
			std::tuple<double, double> GetVelocityAndPosition(void);			//!< 予めセットされたモータ電流と負荷トルクから速度と位置を返す関数(タプル版)
			void SetMotorInertia(const double inertia);				//!< モータ慣性を設定する関数
			void Reset(void);										//!< シミュレータをリセットする関数
			
		private:
			MotorSimulator(const MotorSimulator&) = delete;					//!< コピーコンストラクタ使用禁止
			const MotorSimulator& operator=(const MotorSimulator&) = delete;//!< 代入演算子使用禁止
			
			double Kt;	//!< [Nm/A]	トルク定数
			double Jm;	//!< [kgm^2]	モータ側慣性モーメント
			double Dm;	//!< [Nm/(rad/s)]	モータ側粘性
			double Ts;	//!< [s]	サンプリング時間
			double iq;	//!< [A]	q軸電流
			double taul;//!< [Nm]	負荷トルク
			Matrix<2,2> A;	//!< 連続系 A行列
			Matrix<2,2> B;	//!< 連続系 B行列
			Matrix<1,2> u;	//!< 入力ベクトル
			Matrix<1,2> y;	//!< 出力ベクトル
			StateSpaceSystem<2,2,2> PlantSys;	//!< モータの状態空間モデル
			
			void SetStateSpaceModel(void);	//!< 状態空間モデルを設定する関数
	};
}

#endif

