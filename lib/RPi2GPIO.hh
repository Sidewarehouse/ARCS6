//! @file RPi2GPIO.hh
//! @brief Raspberry Pi 2 GPIOコントローラ
//!
//! ラズベリーパイ2の汎用入出力の制御を行うクラス
//!
//! @date 2019/08/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef RPI2GPIO
#define RPI2GPIO

#include <stdint.h>

namespace ARCS {	// ARCS名前空間
	//! @brief Raspberry Pi 2 GPIOコントローラ
	class RPi2GPIO {
		public:
			//! @brief 入出力モード設定用
			enum IOmode {
				IN,	//!< 入力モード
				OUT	//!< 出力モード
			};
			RPi2GPIO();				//!< コンストラクタ
			~RPi2GPIO();			//!< デストラクタ
			RPi2GPIO(RPi2GPIO&& r);	//!< ムーブコンストラクタ
			void SetConfig(int port, IOmode mode);		// 入出力を設定する関数 port : GPIO番号，mode : RPi2GPIO::IN or RPi2GPIO::OUT
			void SetConfig_GPIO9to0(uint32_t bitdata);	// GPIO9～0の設定レジスタを一括設定する関数
			void SetConfig_GPIO19to10(uint32_t bitdata);// GPIO19～10の設定レジスタを一括設定する関数
			void SetConfig_GPIO29to20(uint32_t bitdata);// GPIO29～20の設定レジスタを一括設定する関数
			void SetConfig_GPIO39to30(uint32_t bitdata);// GPIO39～30の設定レジスタを一括設定する関数
			void SetConfig_GPIO49to40(uint32_t bitdata);// GPIO49～40の設定レジスタを一括設定する関数
			void SetConfig_GPIO53to50(uint32_t bitdata);// GPIO53～50の設定レジスタを一括設定する関数
			void SetConfig_AllOutput(void);				// すべてのGPIOピンを出力モードに設定する関数
			void SetBitDataLo(uint32_t bitdata);		// GPIO31～0 からビットデータを一括出力する関数
			void SetBitDataHi(uint32_t bitdata);		// GPIO39～32 からビットデータを一括出力する関数
			void BitSet(unsigned int port);				// 指定したGPIOを 1 にする関数
			void BitClear(unsigned int port);			// 指定したGPIOを 0 にする関数
			uint32_t GetBitDataLo(void) const;			// GPIO31～0 からデータを一括入力する関数
			uint32_t GetBitDataHi(void) const;			// GPIO53～32 からデータを一括入力する関数
			bool BitGet(unsigned int port) const;		// 指定したGPIOからビットの状態を取得する関数
			void SetACTLED(bool onoff);					// ACT LED の点灯を制御する関数 onoff=trueでオン、onff=falseでオフ
			void SetPWRLED(bool onoff);					// PWR LED の点灯を制御する関数 onoff=trueでオン、onff=falseでオフ
			
		private:
			RPi2GPIO(const RPi2GPIO&) = delete;					//!< コピーコンストラクタ使用禁止
			const RPi2GPIO& operator=(const RPi2GPIO&) = delete;//!< 代入演算子使用禁止
			
			static const int IOREG_SIZE = 64;						//!< [bytes] I/Oレジスタのマッピングサイズ
			static const unsigned long IOREG_ADDRESS = 0x3F200000;	//!< I/Oレジスタのベースアドレス
			
			int fd;			//!< /dev/memファイルディスクリプタ
			uint32_t* IOptr;//!< I/Oレジスタのアドレス
	};
}

#endif

