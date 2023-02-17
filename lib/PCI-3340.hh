//! @file PCI-3340.hh
//! @brief PCI-3340入出力クラス
//!
//! Interface社製PCI-3340のための入出力機能を提供します。(仕様変更済み)
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki & Akira YAMAGUCHI
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PCI_3340
#define PCI_3340

namespace ARCS {	// ARCS名前空間
	//! @brief PCI-3340入出力クラス
	class PCI3340 {
		public:
			explicit PCI3340(unsigned int Base0);	//!< コンストラクタ(DAC初期化＆設定)
			PCI3340();								//!< 空コンストラクタ
			~PCI3340();								//!< デストラクタ(DAC終了処理)
			static const unsigned int MAX_CH = 8;	//!< チャネル最大値
			
			void SetVoltage(double Vout[MAX_CH]);	//!< 指定した電圧[V]をDACから出力
			void SetDigitalOut(unsigned short Val);	//!< 汎用デジタル出力ピンChからValを出力
			unsigned short GetDigitalIn(void);		//!< 汎用デジタル入力ピンChの値を取得
			
		private:
			unsigned int BaseAddr0;					//!< 先頭アドレス0
			
			void Zero(void);									//!< DACの出力電圧を 0[V] にする関数
			void Settings(void);								//!< DACの設定を行う関数
			void Output(unsigned short DACdata[MAX_CH]);		//!< DACから指定した電圧を出力する関数
			static unsigned short IIbyteHi(unsigned short in);	//!< 2byteデータの上位1byteを抽出して出力
			static unsigned short IIbyteLo(unsigned short in);	//!< 2byteデータの下位1byteを抽出して出力
			static unsigned short VoltToDacData(double Vdac);	//!< DAC出力電圧[V]からDACの実際の整数値に変換する
	};
}

#endif



