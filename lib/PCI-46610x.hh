//! @file PCI-46610x.hh
//! @brief Interface社製 PCI-466102～466108 および PEX-466102用クラス
//!
//! RS485/422シリアル通信の送受信機能を提供します。
//! PCI-Express版のPEXボードも使用可能です。
//!
//! @date 2021/09/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PCI_46610X
#define PCI_46610X

#include <stdint.h>
#include <string>

namespace ARCS {	// ARCS名前空間
	//! @brief PCI-466102～466108用クラス(RS485/422シリアル通信ボード)
	class PCI46610x {
		public:
			explicit PCI46610x(unsigned long Addr0);	//!< コンストラクタ
			PCI46610x(PCI46610x&& right);	//!< ムーブコンストラクタ
			~PCI46610x();					//!< デストラクタ
			
			//! @brief ボーレートの定義
			enum BaudRate {
				RATE_2Mbps,		//!< 2Mbps
				RATE_921_6kbps,	//!< 921.6kbps
			};
			
			//! @brief 2線式/4線式の定義
			enum Wire {
				WIRE_2,	//!< 2線式
				WIRE_4,	//!< 4線式
			};
			
			//! @brief パリティの定義
			enum Parity {
				PARITY_DISABLE,	//!< パリティ無効
				PARITY_ODD,		//!< 奇数パリティ
				PARITY_EVEN,	//!< 偶数パリティ
			};
			
			//! @brief ストップビットの定義
			enum StopBit {
				STOPBIT_1,	//!< ストップビット1
				STOPBIT_2,	//!< ストップビット2
			};
			
			//! @brief ワード長の定義
			enum WordLength {
				WORDLENG_5,	//!< ワード長 5bit
				WORDLENG_6,	//!< ワード長 6bit
				WORDLENG_7,	//!< ワード長 7bit
				WORDLENG_8,	//!< ワード長 8bit
			};
			
			static const std::string STR_CR;	//!< キャリッジリターン文字のアスキーコードの定義
			static const std::string STR_LF;	//!< ラインフィード文字のアスキーコードの定義
			
			void SetConfig(
				enum BaudRate rate, enum Wire wir, enum Parity par,
				enum StopBit stpbit, enum WordLength wlen
			);	//!< 通信設定の関数
			void SetTransData(uint8_t TxData);			//!< 送信データ1byte分をセットする関数
			uint8_t GetRecvData(void) const;			//!< 受信データ1byte分を取得する関数
			void GetAllRecvData(std::string& RxData) const;	//!< 全受信データを取得する関数
			unsigned int GetReceiveCount(void) const;	//!< 受信カウント数を返す関数
			void ClearFIFOMemory(void);					//!< 送受信FIFOメモリをクリアする関数
			void ClearAllSettings(void);				//!< すべての設定レジスタをクリアする関数
			
		private:
			PCI46610x(const PCI46610x&) = delete;					//!< コピーコンストラクタ使用禁止
			const PCI46610x& operator=(const PCI46610x&) = delete;	//!< 代入演算子使用禁止
			static constexpr int IOREG_SIZE = 64;					//!< [bytes] I/Oレジスタのマッピングサイズ
			static constexpr unsigned int ADDR_HOLDREG = 0x00;		//!< 送受信ホールディングレジスタのオフセットアドレス
			static constexpr unsigned int ADDR_FIFO_CTRL = 0x02;	//!< FIFOコントロールレジスタのオフセットアドレス
			static constexpr unsigned int ADDR_LINE_CTRL = 0x03;	//!< ラインコントロールレジスタのオフセットアドレス
			static constexpr unsigned int ADDR_LINE_STATUS = 0x05;	//!< ラインステータスレジスタのオフセットアドレス
			static constexpr unsigned int ADDR_RXFIFO_CNT_LO = 0x10;//!< 受信FIFOカウンタ下位のオフセットアドレス
			static constexpr unsigned int ADDR_RXFIFO_CNT_HI = 0x11;//!< 受信FIFOカウンタ上位のオフセットアドレス
			static constexpr unsigned int ADDR_BASECLOCK = 0x20;	//!< 基準クロック周波数設定レジスタのオフセットアドレス
			static constexpr unsigned int ADDR_RESET_REG = 0x21;	//!< ソフトウェアリセットレジスタのオフセットアドレス
			static constexpr unsigned int ADDR_PRESCALER_DLL = 0x28;//!< プリスケーラ(分周器設定の下位8bit)のオフセットアドレス
			static constexpr unsigned int ADDR_PRESCALER_DLM = 0x29;//!< プリスケーラ(分周器設定の上位8bit)のオフセットアドレス
			static constexpr unsigned int ADDR_DRIVE_CTRL = 0x30;	//!< ドライブ制御レジスタのオフセットアドレス
			int fd;					//!< /dev/memファイルディスクリプタ
			uint8_t* IOptr;			//!< I/Oレジスタのアドレス
			unsigned long ADDR_BASE;//!< ベースアドレス
	};
}

#endif

