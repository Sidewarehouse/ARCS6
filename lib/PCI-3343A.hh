//! @file PCI-3343A.hh
//! @brief PCI-3343A入出力クラス
//! Interface社製PCI-3343Aのための入出力機能を提供します。
//! 
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#ifndef PCI_3343A
#define PCI_3343A

namespace ARCS {
	//! @brief PCI-3343A入出力クラス
	//! Interface社製PCI-3343Aのための入出力機能を提供します。
	class PCI3343A {
		public:
			explicit PCI3343A(unsigned int Addr);			//!< コンストラクタ(全チャネル使用する版)
			PCI3343A(unsigned int Addr, uint8_t EnableCh);	//!< コンストラクタ(指定チャネルのみ使用する版)
			PCI3343A();										//!< 空コンストラクタ
			~PCI3343A();									//!< デストラクタ
			void SetVoltage(double V1, double V2, double V3, double V4);	//!< 指定した電圧を出力する関数
			//! V1 : CH1の電圧値[V]  V2 : CH2の電圧値[V]  V3 : CH3の電圧値[V]  V4 : CH4の電圧値[V]
		
		private:
			PCI3343A(const PCI3343A&) = delete;					//!< コピーコンストラクタ使用禁止
			const PCI3343A& operator=(const PCI3343A&) = delete;//!< 代入演算子使用禁止
			
			const unsigned int ADDR_BASE;		//!< ベースアドレス
			const unsigned int ADDR_DACDATA_LO;	//!< DACデータ下位アドレス
			const unsigned int ADDR_DACDATA_HI;	//!< DACデータ上位アドレス
			const unsigned int ADDR_CHSET;		//!< 出力チャンネル設定
			const unsigned int ADDR_CONVMODE;	//!< 変換モード
			const unsigned int ADDR_OUTMODE;	//!< 出力モード
			const unsigned int ADDR_DIO;		//!< 汎用入出力
			uint8_t ENA;						//!< チャネルイネーブル信号 0b0000[Ch4][Ch3][Ch2][Ch1] の並びでビットが「1」 のチャネルが有効
			
			void SetAllEnable(bool flag);		// 全チャネル同時出力イネーブル
			void SetOutEnable(bool flag);		// 全チャネル出力有効 flag = true で電圧出力，flag = false でハイインピーダンス
			void ExecOutput(void);				// 全チャネル同時出力実行(電圧更新)
			void SelectCH(unsigned int ch);		// チャネル選択
			void SetDACdata(uint16_t data);		// DACデータをセットする関数
			void SetAllZero(void);				// 全チャネル零電圧出力
			static uint16_t VoltToDacData(double Vout);	// DAC出力電圧[V]からDACの実際の整数値に変換する関数
			static uint8_t Get2byteHi(uint16_t in);		// 2byteデータの上位1byteを抽出して出力
			static uint8_t Get2byteLo(uint16_t in);		// 2byteデータの下位1byteを抽出して出力
	};
}

#endif
