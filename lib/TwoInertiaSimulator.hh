//! @file TwoInertiaSimulator.hh
//! @brief 2慣性共振系シミュレータV3
//!
//! 2慣性共振系を模擬して，負荷側速度，ねじれ角，モータ速度/位置を計算して出力する。
//!
//! @date 2022/04/06
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef TWOINERTIASIMULATOR
#define TWOINERTIASIMULATOR

#include "Matrix.hh"
#include "StateSpaceSystem.hh"
#include "Integrator.hh"

namespace ARCS {	// ARCS名前空間
//! @brief 2慣性共振系シミュレータV3
class TwoInertiaSimulator {
	public:
		TwoInertiaSimulator(void);	//!< 空コンストラクタ
		TwoInertiaSimulator(const struct TwoInertiaParamDef& Params, const double SmplTime);//!< コンストラクタ
		TwoInertiaSimulator(const struct TwoInertiaParams& Params, const double SmplTime);	//!< コンストラクタ
		TwoInertiaSimulator(TwoInertiaSimulator&& r);	//!< ムーブコンストラクタ
		~TwoInertiaSimulator();							//!< デストラクタ
		void SetParameters(const struct TwoInertiaParamDef& Params, const double SmplTime);	//!< 2慣性共振系のパラメータを設定する関数
		void SetParameters(const struct TwoInertiaParams& Params, const double SmplTime);	//!< 2慣性共振系のパラメータを設定する関数
		void SetCurrentAndLoadTorque(const double Current, const double LoadTorque);		//!< 予めモータ電流と負荷トルクを設定する関数
		void GetResponses(double& LoadSpeed, double& TorsionAngle, double& MotorSpeed);		//!< 予め設定された電流と負荷から負荷側速度＆ねじれ角＆モータ側速度を計算する関数
		void GetResponses(
			const double Current, const double LoadTorque, double& LoadSpeed, double& TorsionAngle, double& MotorSpeed
		);	//!< 負荷側速度＆ねじれ角＆モータ側＆を計算する関数
		std::tuple<double,double,double> GetResponses(double Current, double LoadTorque);	//!< 負荷側速度＆ねじれ角＆モータ側＆を計算する関数(タプル版)
		double GetMotorPosition(void);		//!< モータ側位置応答を取得する関数
		double GetLoadPosition(void);		//!< 負荷側位置を取得する関数
		void SetLoadInertia(double Inertia);//!< 負荷側慣性を設定する関数
		void SetSpring(double Stiffness);	//!< ばね定数を設定する関数
		void Reset(void);					//!< シミュレータをリセットする関数
		
	private:
		TwoInertiaSimulator(const TwoInertiaSimulator&) = delete;					//!< コピーコンストラクタ使用禁止
		const TwoInertiaSimulator& operator=(const TwoInertiaSimulator&) = delete;	//!< 代入演算子使用禁止
		
		double Jl;	//!< [kgm^2]	 負荷側慣性
		double Dl;	//!< [Nm/(rad/s)]負荷側粘性
		double Ds;	//!< [Nm/(rad/s)]2慣性間の粘性
		double Ks;	//!< [Nm/rad]	 2慣性間の剛性
		double Jm;	//!< [kgm^2]	 モータ側慣性
		double Dm;	//!< [Nm/(rad/s)]モータ側粘性
		double Rg;	//!< [-]		 減速比
		double Kt;	//!< [Nm/A]		 トルク定数
		double Ts;	//!< [s]		 サンプリング時間
		Matrix<1,2> u;						//!< 入力ベクトル
		Matrix<1,3> y;						//!< 出力ベクトル
		StateSpaceSystem<3,2,3> PlantSys;	//!< 2慣性系システム
		Integrator<> MotorInteg;			//!< モータ側積分器
		Integrator<> LoadInteg;				//!< 負荷側積分器
		
		void SetStateSpaceModel(void);		//!< 状態空間モデルを設定する関数
};
}

#endif



