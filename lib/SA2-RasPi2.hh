//! @file SA2-RasPi2.hh
//! @brief ServoAmplifier2-RaspberryPi2インターフェースクラス
//!
//! ServoAmplifier2を制御するためのラズベリーパイ2の色々な入出力機能を提供します。
//! (即値マジックナンバーが埋め込まれているが，そのうちキレイにする予定)
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SA2RASPI2
#define SA2RASPI2

#include "RPi2GPIO.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief ServoAmplifier2-RaspberryPi2インターフェースクラス
	class SA2RasPi2 {
		public:
			SA2RasPi2();				//!< コンストラクタ
			SA2RasPi2(SA2RasPi2&& r);	//!< ムーブコンストラクタ
			~SA2RasPi2();				//!< デストラクタ
			void StartProcess(void);			//!< 開始処理
			void EndProcess(void);				//!< 終了処理
			void SetVoltage(double Vu, double Vv, double Vw, double Vdc);	//!< PWMの電圧指令を設定する関数
			void SetDutyRatioPWM(double DutyU, double DutyV, double DutyW);	//!< PWMのデューティ比を設定する関数
			void SetDutyRatioPWMU(double Duty);	//!< PWM U相のデューティ比を設定する関数
			void SetDutyRatioPWMV(double Duty);	//!< PWM V相のデューティ比を設定する関数
			void SetDutyRatioPWMW(double Duty);	//!< PWM W相のデューティ比を設定する関数
			void SetSinglePhaseDutyRatio(double Duty);		//!< 単相Hブリッジインバータとして使う場合のデューティ比を設定する関数
			void SetVoltageDAC(double Vdac1, double Vdac2);	//!< D/Aコンバータの電圧を設定する関数
			void SetVoltageDAC1(double Vdac);				//!< D/AコンバータCH1の電圧を設定する関数
			void SetVoltageDAC2(double Vdac);				//!< D/AコンバータCH2の電圧を設定する関数
			void GetPhaseCurrent(double& Iu, double& Iv, double& Iw);	//!< 相電流を取得する関数
			void GetSinglePhaseCurrent(double& I);	//!< 単相電流を取得する関数
			void GetPosition(double& theta_rm, double& theta_re);		//!< モータ機械角と電気角を取得する関数
			uint32_t GetIncEncoder(void);			//!< インクリメンタルエンコーダのバイナリデータを返す関数
			int32_t GetIncMultiTurnEncoder(void);	//!< インクリメンタルエンコーダの多回転バイナリデータを返す関数
			double GetIncPosition(void);			//!< インクリメンタルエンコーダの位置応答を返す関数
			
		private:
			SA2RasPi2(const SA2RasPi2&) = delete;					//!< コピーコンストラクタ使用禁止
			const SA2RasPi2& operator=(const SA2RasPi2&) = delete;	//!< 代入演算子使用禁止
			static const uint32_t ADDR_PWMU = 0b00000000;		//!< PWM U相 のアドレス
			static const uint32_t ADDR_PWMV = 0b00000001;		//!< PWM V相 のアドレス
			static const uint32_t ADDR_PWMW = 0b00000010;		//!< PWM W相 のアドレス
			static const uint32_t ADDR_DAC1 = 0b00000011;		//!< DAC CH1 のアドレス
			static const uint32_t ADDR_DAC2 = 0b00000100;		//!< DAC CH2 のアドレス
			static const uint32_t ADDR_ABSL = 0b00000000;		//!< 絶対値エンコーダ下位 のアドレス
			static const uint32_t ADDR_ABSH = 0b00000001;		//!< 絶対値エンコーダ上位 のアドレス
			static const uint32_t ADDR_IU   = 0b00000010;		//!< 電流センサ U相 のアドレス
			static const uint32_t ADDR_IV   = 0b00000011;		//!< 電流センサ V相 のアドレス
			static const uint32_t ADDR_IW   = 0b00000100;		//!< 電流センサ W相 のアドレス
			static const uint32_t ADDR_INCL = 0b00000101;		//!< インクリメンタルエンコーダ下位 のアドレス
			static const uint32_t ADDR_INCH = 0b00000110;		//!< インクリメンタルエンコーダ上位 のアドレス
			static const uint32_t INC_MASK  = 0x00000FFF;		//!< インクリメンタルエンコーダのマスク(分解能によって変える)
			RPi2GPIO gpio;		//!< ラズパイGPIO
			uint32_t DATA_PWMU;	//!< PWM U相 のデータバッファ
			uint32_t DATA_PWMV;	//!< PWM V相 のデータバッファ
			uint32_t DATA_PWMW;	//!< PWM W相 のデータバッファ
			uint32_t DATA_DAC1;	//!< DAC CH1 のデータバッファ
			uint32_t DATA_DAC2;	//!< DAC CH2 のデータバッファ
			uint32_t DATA_IU;	//!< 電流センサ U相 のデータバッファ
			uint32_t DATA_IV;	//!< 電流センサ V相 のデータバッファ
			uint32_t DATA_IW;	//!< 電流センサ W相 のデータバッファ
			uint32_t DATA_ABS;	//!< 絶対値エンコーダのデータバッファ
			uint32_t DATA_INC;	//!< インクリメンタルエンコーダのデータバッファ
			uint32_t RawCountPrev;	//!< 絶対値多回転カウント値計算用のバッファ
			int32_t UpperCount;		//!< 絶対値多回転カウント値計算用の上位カウンタ
			int32_t AbsMultCount;	//!< 絶対値多回転カウント値
			uint32_t IncRawCountPrev;	//!< インクリメンタル多回転カウント値計算用のバッファ
			int32_t IncUpperCount;		//!< インクリメンタル多回転カウント値計算用の上位カウンタ
			void WaitPwmSyncHigh(void);							//!< PWM_SYNC が High になるまでブロッキング待機する関数
			void WaitPwmSyncLow(void);							//!< PWM_SYNC が Low になるまでブロッキング待機する関数
			void WaitEndOfConversionIuvw(void);					//!< 電流センサのAD変換が完了するまでブロッキング待機する関数
			void WaitEndOfConversionAbs(void);					//!< 絶対値エンコーダのSP変換が完了するまでブロッキング待機する関数
			void WaitEndOfConversionInc(void);					//!< インクリメンタルエンコーダの変換処理が完了するまでブロッキング待機する関数
			void SendConversionStart(void);						//!< CONVST信号を出力する関数
			void SendWriteSelect(void);							//!< WR_SEL信号を出力する関数
			void SetReadSelect(void);							//!< RD_SEL信号を出力する関数
			void ClearReadSelect(void);							//!< RD_SEL信号を停止する関数
			void SetBusModeIn(void);							//!< データバスを入力モードに設定する関数
			void SetBusModeOut(void);							//!< データバスを出力モードに設定する関数
			void WaitUsingNop(const unsigned int N) const;		//!< NOP命令を使って待機する関数
			void WriteDataBus(uint32_t Address, uint32_t Data);	//!< 指定したアドレスでデータバスに書き込む関数
			uint32_t ReadDataBus(uint32_t Address);				//!< 指定したアドレスでデータバスから読み込む関数
			int32_t GetMultiTurnCount(uint32_t RawCount);		//!< 絶対値多回転カウント値を取得する関数
			int32_t GetIncMultiTurnCount(uint32_t RawCount);	//!< インクリメンタル多回転カウント値を取得する関数
			static uint32_t VoltToDacData(double Vdac);			//!< 電圧[V]からDACバイナリデータへの換算関数
			static uint32_t DutyToPwmData(double Duty);			//!< デューティ比からPWMバイナリデータへの換算関数
			static double AdcDataToAmpere(uint32_t AdcData);	//!< ADCバイナリデータから電流[A]への換算関数
			static double AbsDataToRadian(int32_t AbsMultData);	//!< エンコーダデータから位置[rad]への換算関数
			static double ConvElectAngle(uint32_t AbsData);		//!< モータ電気角 [rad] へ換算する関数
	};
}

#endif

