//! @file InterfaceFunctions.hh
//! @brief インターフェースクラス
//! @date 2024/06/25
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef INTERFACEFUNCTIONS
#define INTERFACEFUNCTIONS

// 基本のインクルードファイル
#include <array>
#include "ConstParams.hh"
#include "ARCSeventlog.hh"
#include "ARCSassert.hh"
#include "ARCSprint.hh"
#include "EquipParams.hh"

// 追加のARCSライブラリをここに記述
#include "ArcsMatrix.hh"
#include "Limiter.hh"

namespace ARCS {	// ARCS名前空間
//! @brief インターフェースクラス
//! 「電流指令,位置,トルク,…等々」と「DAC,エンコーダカウンタ,ADC,…等々」との対応を指定します。
class InterfaceFunctions {
	public:
		// ここにインターフェース関連の定数を記述する(記述例はsampleを参照)
		
		// ここにD/A，A/D，エンコーダIFボードクラス等々の宣言を記述する(記述例はsampleを参照)
		
		//! @brief コンストラクタ
		InterfaceFunctions()
			// ここにD/A，A/D，エンコーダIFボードクラス等々の初期化子リストを記述する(記述例はsampleを参照)
			//:
		{
			PassedLog();
		}

		//! @brief デストラクタ
		~InterfaceFunctions(){
			SetZeroCurrent();	// 念のためのゼロ電流指令
			PassedLog();
		}

		//! @brief サーボON信号を送出する関数
		void ServoON(void){
			// ここにサーボアンプへのサーボON信号の送出シーケンスを記述する
			
		}

		//! @brief サーボOFF信号を送出する関数
		void ServoOFF(void){
			// ここにサーボアンプへのサーボOFF信号の送出シーケンスを記述する
			
		}
		
		//! @brief 電流指令をゼロに設定する関数
		void SetZeroCurrent(void){
			// ここにゼロ電流指令とサーボアンプの関係を列記する
			
		}
		
		//! @brief 位置ベクトルを取得する関数
		//! @param[out]	Position	位置ベクトル [rad]
		void GetPosition(ArcsMat<EquipParams::ACTUATOR_NUM, 1>& Position){
			// ここにエンコーダとPositionベクトルとの関係を列記する
			
		}
		
		//! @brief 位置と速度ベクトルを取得する関数
		//! @param[out]	Position	位置ベクトル [rad]
		//! @param[out]	Velocity	速度ベクトル [rad/s]
		void GetPositionAndVelocity(
			ArcsMat<EquipParams::ACTUATOR_NUM, 1>& Position,
			ArcsMat<EquipParams::ACTUATOR_NUM, 1>& Velocity
		){
			// ここにエンコーダ，速度演算結果とPosition，Velocityベクトルとの関係を列記する
			
		}
		
		//! @brief モータ電気角と機械角ベクトルを取得する関数
		//! @param[out]	ElectAngle	電気角ベクトル [rad]
		//! @param[out]	MechaAngle	機械角ベクトル [rad]
		void GetElectricAndMechanicalAngle(
			ArcsMat<EquipParams::ACTUATOR_NUM, 1>& ElectAngle,
			ArcsMat<EquipParams::ACTUATOR_NUM, 1>& MechaAngle
		){
			// ここにモータ電気角，機械角とElecAngle，MechaAngleベクトルとの関係を列記する
			
		}
		
		//! @brief トルクベクトルを取得する関数
		//! @param[out]	Torque	トルクベクトル [Nm]
		void GetTorque(ArcsMat<EquipParams::ACTUATOR_NUM, 1>& Torque){
			// ここにトルクセンサとTorqueベクトルとの関係を列記する
			
		}
		
		//! @brief 加速度応答を取得する関数
		//! @param[out]	Acceleration	加速度ベクトル [rad/s^2]
		void GetAcceleration(ArcsMat<3,1>& Acceleration){
			// ここに加速度センサとAccelerationベクトルとの関係を列記する
			
		}
		
		//! @brief 電流ベクトルを取得する関数
		//! @param[out]	Current	電流ベクトル [A]
		void GetCurrent(ArcsMat<EquipParams::ACTUATOR_NUM, 1>& Current){
			// ここに電流センサとCurrentベクトルとの関係を列記する
			
		}
		
		//! @brief 電流指令ベクトルを設定する関数
		//! @param[in]	CurrentRef	電流指令ベクトル [A]
		void SetCurrent(const ArcsMat<EquipParams::ACTUATOR_NUM, 1>& CurrentRef){
			// ここにCurrentRefベクトルとサーボアンプの関係を列記する
			
		}
		
		//! @brief トルク指令ベクトルを設定する関数
		//! @param[in]	TorqueRef	トルク指令ベクトル [Nm]
		void SetTorque(const ArcsMat<EquipParams::ACTUATOR_NUM, 1>& TorqueRef){
			// ここにTorqueRefベクトルとサーボアンプの関係を列記する
			
		}
		
		//! @brief 6軸力覚センサ応答を取得する関数
		//! @param[out]	ForceTorque [ XYZ軸並進力[N], rpw軸トルク[Nm] ]^T のベクトル
		void Get6axisForce(ArcsMat<6,1>& ForceTorque){
			// ここに6軸力覚センサと力/トルクベクトルとの関係を列記する
			
		}
		
		//! @brief 安全装置への信号出力を設定する関数
		//! @param[in]	Signal	安全装置へのディジタル信号
		void SetSafetySignal(const uint8_t& Signal){
			// ここに安全信号とDIOポートとの関係を列記する
			
		}
		
		//! @brief Z相クリアに関する設定をする関数
		//! @param[in]	ClearEnable	true = Z相が来たらクリア，false = クリアしない
		void SetZpulseClear(const bool ClearEnable){
			// インクリメンタルエンコーダのZ(I,C)相クリアの設定が必要な場合に記述する
			
		}
		
	private:
		InterfaceFunctions(const InterfaceFunctions&) = delete;					//!< コピーコンストラクタ使用禁止
		const InterfaceFunctions& operator=(const InterfaceFunctions&) = delete;//!< 代入演算子使用禁止
		
		// ここにセンサ取得値とSI単位系の間の換算に関する関数を記述(記述例はsampleを参照)
		
		//! @brief モータ機械角 [rad] へ換算する関数
		//! @brief	count	エンコーダカウント値
		//! @return	機械角 [rad]
		static double ConvMotorAngle(const long count){
			return 0;	//ENC_TO_RADIAN*(double)count;
		}
		
		//! @brief モータ電気角 [rad] へ換算する関数 (-2π～+2πの循環値域制限あり)
		//! @brief	count	エンコーダカウント値
		//! @return	電気角 [rad]
		static double ConvElectAngle(const long count){
			return 0;	//ENC_TO_RADIAN*(double)(ENC_POLEPARE*( count % (ENC_MAX_COUNT/ENC_POLEPARE) ));
		}
};
}

#endif

