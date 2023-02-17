//! @file PCI-3133.hh
//! @brief PCI-3133入出力クラス
//!
//! Interface社製PCI-3133のための入出力機能を提供します。(1チャネルのみ実装済み)
//!
//! @date 2019/02/26
//! @author HIDETAKA MORIMITSU & Yuki YOKOKURA
//
// Copyright (C) 2011-2019 HIDETAKA MORIMITSU & Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PCI_3133
#define PCI_3133

namespace ARCS {	// ARCS名前空間
	//! @brief PCI-3133入出力クラス
	class PCI3133 {
		public:
			explicit PCI3133(const unsigned int Base);	//!< コンストラクタ(ADC初期化＆設定)
			PCI3133();									//!< 空コンストラクタ(ADC初期化＆設定)
			~PCI3133();									//!< デストラクタ(ADC終了処理)
			double GetVoltage(void);					//!< 電圧値[V]を取得
			void Settings(void);						//!< A/Dコンバータの設定を行う関数
			
		private:
			unsigned int BaseAddr;							//!< 先頭アドレス
			unsigned short Input(void);						//!< A/Dコンバータから電圧データを読み込む関数
			static double AdcDataToVolt(unsigned short ADCdata);	//!< 電圧データから電圧値[V]に変換する関数
			static unsigned long IIbyteCat(unsigned short High, unsigned short Low);	//!< 上位、下位に分かれている各々1バイトのデータを、2バイトのデータに結合する
	};
}

#endif



