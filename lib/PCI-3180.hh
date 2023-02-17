//! @file PCI-3180.hh
//! @brief PCI-3180入出力クラス
//! Interface社製PCI-3180のための入出力機能を提供します。
//! @date 2020/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#ifndef PCI_3180
#define PCI_3180

namespace ARCS {
	//! @brief PCI-3180入出力クラス
	//! Interface社製PCI-3180のための入出力機能を提供します。
	class PCI3180 {
		public:
			//! @brief 入力電圧レンジ定義
			enum RangeMode {
				RANGE_B_10V,	//!< バイポーラ±10V入力モード
				RANGE_B_5V,		//!< バイポーラ ±5V入力モード
				RANGE_U_10V		//!< ユニポーラ  10V入力モード
			};
			
			PCI3180(unsigned int Addr, RangeMode Range);	//!< コンストラクタ(全チャネルを使用する版)
			PCI3180(unsigned int Addr, RangeMode Range, uint8_t EnableCh);	//!< コンストラクタ(指定したチャネルのみを使用する版)
			PCI3180();				//!< 空コンストラクタ
			~PCI3180();				//!< デストラクタ
			
			void ConvStart(void);	//!< AD変換開始
			bool GetBusy(void);		//!< 変換待機信号の取得 true = 変換中, false = 変換終了
			void WaitBusy(void);	//!< AD変換が完了するまでポーリング待機(ブロッキング動作関数)
			void GetVoltage(double& V1, double& V2, double& V3, double& V4);	//!< 電圧を取得する関数
			void SetDigitalOutput(const uint8_t& Data);								//!< ディジタル出力ポートを設定する関数
		
		private:
			PCI3180(const PCI3180&) = delete;					//!< コピーコンストラクタ使用禁止
			const PCI3180& operator=(const PCI3180&) = delete;	//!< 代入演算子使用禁止
			
			const unsigned int ADDR_BASE;		//!< ベースアドレス
			const unsigned int ADDR_ADCDATA_LO;	//!< ADCデータ下位アドレス
			const unsigned int ADDR_ADCDATA_HI;	//!< ADCデータ上位アドレス
			const unsigned int ADDR_CHSET;		//!< 入力チャネル設定
			const unsigned int ADDR_BUSY;		//!< 変換待機信号
			const unsigned int ADDR_CONVST;		//!< 変換開始信号
			const unsigned int ADDR_RANGE;		//!< 入力レンジ設定
			const unsigned int ADDR_RANGE_CH;	//!< 入力レンジ設定チャンネル選択
			const unsigned int ADDR_DIO;		//!< デジタル出力(2ch)
			
			double A;	//!< 電圧換算の傾き
			double B;	//!< 電圧換算の切片
			uint8_t ENA;//!< チャネルイネーブル信号 0b0000[Ch4][Ch3][Ch2][Ch1] の並びでビットが「1」 のチャネルが有効
			
			void SelectCH(unsigned int ch);		//!< チャネル選択
			uint16_t GetADCdata(void);			//!< ADCデータを取得する関数
			double AdcDataToVolt(uint16_t data) const;//!< ADCデータから電圧値 [V] に変換する関数
			void SelectRangeCH(unsigned int ch);//!< 入力レンジ設定チャンネル選択
			void SetInputRange(RangeMode Range);//!< 入力レンジの設定
			void CalcVoltConv(RangeMode Range);	//!< 電圧換算の傾きと切片の計算
			static uint16_t Combine2byte(uint8_t High, uint8_t Low);	//!< 上位、下位に分かれている各々1バイトのデータを、2バイトのデータに結合する
	};
}

#endif

