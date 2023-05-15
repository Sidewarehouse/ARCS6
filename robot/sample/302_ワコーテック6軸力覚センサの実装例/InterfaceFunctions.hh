//! @file InterfaceFunctions.hh
//! @brief インターフェースクラス
//! @date 2021/03/31
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#ifndef INTERFACEFUNCTIONS
#define INTERFACEFUNCTIONS

// 基本のインクルードファイル
#include <array>
#include "ConstParams.hh"
#include "ARCSeventlog.hh"
#include "ARCSassert.hh"
#include "ARCSprint.hh"

// 追加のARCSライブラリをここに記述
#include "Limiter.hh"
#include "PCI-46610x.hh"
#include "WEF-6A.hh"

namespace ARCS {	// ARCS名前空間
//! @brief インターフェースクラス
//! 「電流指令,位置,トルク,…等々」と「DAC,エンコーダカウンタ,ADC,…等々」との対応を指定します。
class InterfaceFunctions {
	public:
		// ここにインターフェース関連の定数を記述する(記述例はsampleを参照)
		unsigned long ADDR_RS422CH1 = 0xd0401000;
		
		// ここにD/A，A/D，エンコーダIFボードクラス等々の宣言を記述する(記述例はsampleを参照)
		std::unique_ptr<PCI46610x> RS422CH1;// RS422シリアル通信ボード
		std::unique_ptr<WEF6A> ForceSensor;	// 6軸力覚センサ
		
		//! @brief コンストラクタ
		InterfaceFunctions()
			// ここにD/A，A/D，エンコーダIFボードクラス等々の初期化子リストを記述する(記述例はsampleを参照)
			: RS422CH1(nullptr), ForceSensor(nullptr)
		{
			PassedLog();
			RS422CH1 = std::make_unique<PCI46610x>(ADDR_RS422CH1);	// RS422シリアル通信ボード
			RS422CH1->SetConfig(
				PCI46610x::RATE_921_6kbps,	// ボーレート 921.6kbps
				PCI46610x::WIRE_4,			// 4線式
				PCI46610x::PARITY_DISABLE,	// パリティなし
				PCI46610x::STOPBIT_1,		// ストップビット1
				PCI46610x::WORDLENG_8		// ワード長8bit
			);
			ForceSensor = std::make_unique<WEF6A>(std::move(RS422CH1));
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
		
		//! @brief 位置応答を取得する関数
		//! @param[out]	PositionRes	位置応答 [rad]
		void GetPosition(std::array<double, ConstParams::ACTUATOR_NUM>& PositionRes){
			// ここにエンコーダとPosition配列との関係を列記する
			
		}
		
		//! @brief 位置応答と速度応答を取得する関数
		//! @param[out]	PositionRes	位置応答 [rad]
		//! @param[out]	VelocityRes	速度応答 [rad/s]
		void GetPositionAndVelocity(
			std::array<double, ConstParams::ACTUATOR_NUM>& PositionRes,
			std::array<double, ConstParams::ACTUATOR_NUM>& VelocityRes
		){
			// ここにエンコーダ，速度演算結果とPositionRes配列，VelocityRes配列との関係を列記する
			
		}
		
		//! @brief モータ電気角と機械角を取得する関数
		//! @param[out]	ElectAngle	電気角 [rad]
		//! @param[out]	MechaAngle	機械角 [rad]
		void GetElectricAndMechanicalAngle(
			std::array<double, ConstParams::ACTUATOR_NUM>& ElectAngle,
			std::array<double, ConstParams::ACTUATOR_NUM>& MechaAngle
		){
			// ここにモータ電気角，機械角とElePosition配列，MecPosition配列との関係を列記する
			
		}
		
		//! @brief トルク応答を取得する関数
		//! @param[out]	Torque	トルク応答 [Nm]
		void GetTorque(std::array<double, ConstParams::ACTUATOR_NUM>& Torque){
			// ここにトルクセンサとTorque配列との関係を列記する
			
		}
		
		//! @brief 加速度応答を取得する関数
		//! @param[out]	Acceleration	加速度応答 [rad/s^2]
		void GetAcceleration(std::array<double, ConstParams::ACTUATOR_NUM>& Acceleration){
			// ここに加速度センサとAcceleration配列との関係を列記する
			
		}
		
		//! @brief トルク応答と加速度応答を取得する関数
		//! @param[out]	Torque	トルク応答 [Nm]
		//! @param[out]	Acceleration	加速度応答 [rad/s^2]
		void GetTorqueAndAcceleration(
			std::array<double, ConstParams::ACTUATOR_NUM>& Torque,
			std::array<double, ConstParams::ACTUATOR_NUM>& Acceleration
		){
			// ここにトルクセンサとTorque配列との関係，加速度センサとAcceleration配列との関係を列記する
			
		}
		
		//! @brief 電流応答を取得する関数
		//! @param[out]	Current	電流応答 [A]
		void GetCurrent(std::array<double, ConstParams::ACTUATOR_NUM>& Current){
			// ここに電流センサとCurrent配列との関係を列記する
			
		}
		
		//! @brief 電流指令を設定する関数
		//! @param[in]	Current	電流指令 [A]
		void SetCurrent(const std::array<double, ConstParams::ACTUATOR_NUM>& Current){
			// ここにCurrent配列とサーボアンプの関係を列記する
			
		}
		
		//! @brief トルク指令を設定する関数
		//! @param[in]	Torque	トルク指令 [Nm]
		void SetTorque(const std::array<double, ConstParams::ACTUATOR_NUM>& Torque){
			// ここにTorque配列とサーボアンプの関係を列記する
			
		}
		
		//! @brief 6軸力覚センサ応答を取得する関数
		//! @param[out]	Fx-Fy 各軸の並進力 [N]
		//! @param[out]	Mx-My 各軸のトルク [Nm]
		void Get6axisForce(double& Fx, double& Fy, double& Fz, double& Mx, double& My, double& Mz){
			// ここに6軸力覚センサと各変数との関係を列記する
			ForceSensor->SendForceRequest();	// リクエスト送信
			ForceSensor->WaitForceData();		// データ取得待機(ブロッキング動作)
			ForceSensor->Get6axisForce(Fx, Fy, Fz, Mx, My, Mz);	// [N],[Nm] 力覚センサ値の取得
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

