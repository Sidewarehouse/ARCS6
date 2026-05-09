//! @file MotorSimulator.hh
//! @brief モータシミュレータ
//!
//! モータを模擬する
//!
//! @date 2026/05/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2025 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef MOTORSIMULATOR
#define MOTORSIMULATOR

#include "ArcsMatrix.hh"
#include "ArcsControl.hh"
#include "StateSpaceSystem.hh"

namespace ARCS {	// ARCS名前空間
//! @brief モータシミュレータ
class MotorSimulator {
	public:
		//! @brief 空コンストラクタ
		MotorSimulator(void)
			: Kt(0), Jm(0), Dm(0), Ts(0), iq(0), taul(0),
			A(), B(), u(), y(), Plant()
		{
			PassedLog();
		}

		//! @brief コンストラクタ
		//! @param[in]	TrqConst	[Nm/A] トルク定数
		//! @param[in]	MotorInert	[kgm^2]モータ慣性
		//! @param[in]	MotorFric	[Nm/(rad/s)]モータ粘性
		//! @param[in]	SmplTime	[s] サンプリング時間
		MotorSimulator(const double TrqConst, const double MotorInert, const double MotorFric, const double SmplTime)
			: Kt(TrqConst), Jm(MotorInert), Dm(MotorFric), Ts(SmplTime), iq(0), taul(0),
			A(), B(), u(), y(), Plant()
		{
			SetStateSpaceModel();	// 状態空間モデルにセット
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		MotorSimulator(MotorSimulator&& r)
			: Kt(r.Kt), Jm(r.Jm), Dm(r.Dm), Ts(r.Ts), iq(r.iq), taul(r.taul),
			A(r.A), B(r.B), u(r.u), y(r.y), Plant(std::move(r.Plant))
		{
			
		}

		//! @brief デストラクタ
		~MotorSimulator(){
			PassedLog();
		}

		//! @brief パラメータを設定する関数
		//! @param[in]	TrqConst	[Nm/A] トルク定数
		//! @param[in]	MotorInert	[kgm^2]モータ慣性
		//! @param[in]	MotorFric	[Nm/(rad/s)]モータ粘性
		//! @param[in]	SmplTime	[s] サンプリング時間
		void SetParameters(const double TrqConst, const double MotorInert, const double MotorFric, const double SmplTime){
			Kt = TrqConst;
			Jm = MotorInert;
			Dm = MotorFric;
			Ts = SmplTime;
			SetStateSpaceModel();	// 状態空間モデルにセット
		}

		//! @brief モータ電流と負荷トルクを設定する関数
		//! @param[in]	current	[A] 電流
		//! @param[in]	loadtorque	[Nm] 負荷トルク
		void SetCurrentAndLoadTorque(const double current, const double loadtorque){
			// 入力ベクトルの設定
			u.Set(
				current,
				loadtorque
			);
			Plant.SetInput(u);
		}

		//! @brief モータ速度と位置を返す関数(引数渡し版)
		//! @param[out]	velocity	モータ速度 [rad/s]
		//! @param[out]	position	位置 [rad]
		void GetSpeedAndPosition(double& speed, double& position){
			// 応答計算
			Plant.Update();
			
			// 状態ベクトルから抽出
			std::tie(speed, position) = Plant.GetNextOutput2();
		}

		//! @brief モータ速度と位置を返す関数(タプル返し版)
		//! @return	モータ速度 [rad/s], 位置 [rad]
		std::tuple<double,double> GetSpeedAndPosition(void){
			double speed, position;
			GetSpeedAndPosition(speed, position);	// 応答計算
			
			// 速度と位置をタプルで返す
			return {speed, position};
		}

		//! @brief モータ慣性を設定する関数
		//! @param[in]	inertia	モータ慣性 [kgm^2]
		void SetMotorInertia(const double inertia){
			Jm = inertia;
			SetStateSpaceModel();	// 状態空間モデルにセット
		}

		//! @brief シミュレータをリセットする関数
		void Reset(void){
			Plant.ClearStateVector();
		}

	private:
		MotorSimulator(const MotorSimulator&) = delete;					//!< コピーコンストラクタ使用禁止
		const MotorSimulator& operator=(const MotorSimulator&) = delete;//!< 代入演算子使用禁止
		
		double Kt;	//!< [Nm/A]	トルク定数
		double Jm;	//!< [kgm^2]	モータ側慣性モーメント
		double Dm;	//!< [Nm/(rad/s)]	モータ側粘性
		double Ts;	//!< [s]	サンプリング時間
		double iq;	//!< [A]	q軸電流
		double taul;//!< [Nm]	負荷トルク
		ArcsMat<2,2> A;	//!< 連続系 A行列
		ArcsMat<2,2> B;	//!< 連続系 B行列
		ArcsMat<2,1> u;	//!< 入力ベクトル
		ArcsMat<2,1> y;	//!< 出力ベクトル
		ArcsControl::StateSpace<2,2,2> Plant;	//!< モータの状態空間モデル
		
		//! @brief 状態空間モデルを設定する関数
		void SetStateSpaceModel(void){
			// 連続系A行列の設定
			A.Set(
				-Dm/Jm,  0,
					1,  0
			);
			
			// 連続系B行列の設定
			B.Set(
				Kt/Jm, -1.0/Jm,
					0,       0
			);
			
			auto C = ArcsMat<2,2>::eye();		// C行列の設定
			Plant.SetSystem(A, B, C, Ts);
		}
};
}

#endif

