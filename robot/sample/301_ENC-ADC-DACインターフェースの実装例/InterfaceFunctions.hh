//! @file InterfaceFunctions.hh
//! @brief インターフェースクラス
//! @date 2020/04/08
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
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
#include "PCI-6205C.hh"
#include "PCI-3180.hh"
#include "PCI-3343A.hh"

namespace ARCS {	// ARCS名前空間
//! @brief インターフェースクラス
//! 「電流指令,位置,トルク,…等々」と「DAC,エンコーダカウンタ,ADC,…等々」との対応を指定します。
class InterfaceFunctions {
	public:
		// ここにインターフェース関連の定数を記述する(記述例はsampleを参照)
		// エンコーダ関連の定数
		static constexpr unsigned int PCI6205C_ADDR1 = 0xE880;	// ENC PCI-6205C のベースアドレス設定 ("lspci -vv" で調べて入力)
		static constexpr unsigned int PCI6205C_ADDR2 = 0xE800;
		static constexpr unsigned int PCI6205C_ADDR3 = 0xE480;
		static constexpr unsigned int PCI6205C_ADDR4 = 0xE400;
		static constexpr bool ENC_MULT4 = true;					// 4逓倍設定
		static constexpr long ENC_MAX_COUNT = 8000;				// [pulse] エンコーダ最大値(逓倍後)
		static constexpr int ENC_POLEPARE = 5;					// [-] 極対数
		static constexpr double ENC_TO_RADIAN = 2.0*M_PI/((double)ENC_MAX_COUNT);// [rad/pulse] ラジアンへの換算値
		
		// トルクセンサ＆加速度センサ関連の定数
		static constexpr unsigned int PCI3180_ADDR = 0xE000;	// ADC PCI-3180 のベースアドレス設定 ("lspci -vv" で調べて入力)
		static constexpr unsigned int PCI3180_ENA  = 0b00000111;// ADC イネーブルにするチャネルの指定
		static constexpr double TRQSEN_TO_NM = -6;				// [Nm/V] トルクへの換算値
		static constexpr double ACCSEN_TO_MPS2 = 49.02;			// [(m/s^2)/V] 接線加速度への換算値
		static constexpr double ACCSEN_RADIUS = 0.2;			// [m] 加速度センサ設置半径
		
		// サーボアンプ関連の定数
		static constexpr unsigned int PCI3343A_ADDR= 0xEC00;		// DAC PCI-3343A のベースアドレス設定 ("lspci -vv" で調べて入力)
		static constexpr unsigned int PCI3343A_ENA = 0b00001111;	// DAC イネーブルにするチャネルの指定
		static constexpr double SRV_A_TO_V = 5.0/1.3*sqrt(2.0/3.0);	// [V/A] 電流電圧換算ゲイン
		
		// ここにD/A，A/D，エンコーダIFボードクラス等々の宣言を記述する(記述例はsampleを参照)
		PCI6205C ENC;
		PCI3180  ADC;
		PCI3343A DAC;
		
		//! @brief コンストラクタ
		InterfaceFunctions()
			// ここにD/A，A/D，エンコーダIFボードクラス等々の初期化子リストを記述する(記述例はsampleを参照)
			: ENC(PCI6205C_ADDR1, PCI6205C_ADDR2, PCI6205C_ADDR3, PCI6205C_ADDR4, ConstParams::ACTUATOR_NUM, ENC_MULT4),
			  ADC(PCI3180_ADDR, PCI3180::RANGE_B_5V, PCI3180_ENA),
			  DAC(PCI3343A_ADDR, PCI3343A_ENA)
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
			DAC.SetVoltage(0,0,0,0);
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
			std::array<long, PCI6205C::MAX_CH> count = {0};
			ENC.GetCount(count);	// エンコーダカウント数を取得
			// モータ電気角の取得
			ElectAngle[0] = ConvElectAngle( count[0]);	// [rad] 1軸目
			ElectAngle[1] = ConvElectAngle(-count[1]);	// [rad] 2軸目 (エンコーダが逆に取り付けられているので注意)
			// モータ機械角の取得
			MechaAngle[0] = ConvMotorAngle( count[0]);	// [rad] 1軸目
			MechaAngle[1] = ConvMotorAngle(-count[1]);	// [rad] 2軸目 (エンコーダが逆に取り付けられているので注意)
		}
		
		//! @brief トルク応答と加速度応答を取得する関数
		//! @param[out]	Torque	トルク応答 [Nm]
		//! @param[out]	Acceleration	加速度応答 [rad/s^2]
		void GetTorqueAndAcceleration(
			std::array<double, ConstParams::ACTUATOR_NUM>& Torque,
			std::array<double, ConstParams::ACTUATOR_NUM>& Acceleration
		){
			// ここにトルクセンサとTorque配列との関係，加速度センサとAcceleration配列との関係を列記する
			double V1, V2, V3, V4;				// [V] 入力電圧
			ADC.ConvStart();					// AD変換開始
			ADC.WaitBusy();						// AD変換が完了するまで待機(ブロッキング動作)
			ADC.GetVoltage(V1, V2, V3, V4);		// [V] 電圧値の取得
			Torque[0] = TRQSEN_TO_NM*V1;		// [Nm] 電圧からモータ1のトルク応答へ換算
			Torque[1] = TRQSEN_TO_NM*V2;		// [Nm] 電圧からモータ2のトルク応答へ換算
			Acceleration[0] = ACCSEN_TO_MPS2*V3/ACCSEN_RADIUS;	// [rad/s^2] 電圧からモータ1の負荷側加速度応答へ換算
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
			double V1, V2, V3, V4;			// [V] 出力電圧
			V1 = -SRV_A_TO_V*Limiter(Current[0], ConstParams::ACT_MAX_CURRENT[0]);	// [V] U相電流指令 モータ1 電流リミッタをかけてから電圧信号へ換算
			V2 = -SRV_A_TO_V*Limiter(Current[1], ConstParams::ACT_MAX_CURRENT[0]);	// [V] V相電流指令 モータ1
			V3 = -SRV_A_TO_V*Limiter(Current[2], ConstParams::ACT_MAX_CURRENT[1]);	// [V] U相電流指令 モータ2
			V4 = -SRV_A_TO_V*Limiter(Current[3], ConstParams::ACT_MAX_CURRENT[1]);	// [V] V相電流指令 モータ2
			DAC.SetVoltage(V1, V2, V3, V4);	// DA変換電圧更新
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
			ENC.ZpulseClear(ClearEnable);
		}
		
	private:
		InterfaceFunctions(const InterfaceFunctions&) = delete;					//!< コピーコンストラクタ使用禁止
		const InterfaceFunctions& operator=(const InterfaceFunctions&) = delete;//!< 代入演算子使用禁止
		
		// ここにセンサ取得値とSI単位系の間の換算に関する関数を記述(記述例はsampleを参照)
		
		//! @brief モータ機械角 [rad] へ換算する関数
		//! @brief	count	エンコーダカウント値
		//! @return	機械角 [rad]
		static double ConvMotorAngle(const long count){
			return ENC_TO_RADIAN*(double)count;
		}
		
		//! @brief モータ電気角 [rad] へ換算する関数 (-2π～+2πの循環値域制限あり)
		//! @brief	count	エンコーダカウント値
		//! @return	電気角 [rad]
		static double ConvElectAngle(const long count){
			return ENC_TO_RADIAN*(double)(ENC_POLEPARE*( count % (ENC_MAX_COUNT/ENC_POLEPARE) ));
		}
};
}

#endif

