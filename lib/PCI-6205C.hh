//! @file PCI-6205C.hh
//! @brief PCI-6205C入出力クラス
//!
//! Interface社製PCI-6205Cのための入出力機能を提供します。
//!
//! @date 2020/04/08
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef PCI_6205C
#define PCI_6205C

#include <sys/io.h>
#include <unistd.h>
#include <stdint.h>
#include <algorithm>
#include <array>

namespace ARCS {	// ARCS名前空間
	//! @brief PCI-6205C入出力クラス
	class PCI6205C {
		public:
			PCI6205C(const unsigned int Addr0, const unsigned int Addr1, const unsigned int Addr2, const unsigned int Addr3,
						unsigned int NumOfChannel, bool Mult4Mode);	//!< コンストラクタ(ENC初期化＆設定) Addr0-3 : アドレス, NumOfChannel : 仕様チャネル数, Mult4Mode : 4逓倍するか否か
			PCI6205C();								//!< 空コンストラクタ
			~PCI6205C();							//!< デストラクタ(ENC終了処理)
			static const unsigned int MAX_CH = 8;	//!< チャネル最大値
			
			void ZpulseClear(bool flag);			//!< Z相クリア設定
			void GetCount(long count[MAX_CH]);		//!< エンコーダカウンタからカウント値を読み込む関数(生配列版)
			void GetCount(std::array<long, MAX_CH>& count);	//!< エンコーダカウンタからカウント値を読み込む関数(std::array版)
			void ClearAllCounter(void);				//!< エンコーダカウンタの値を零(0x000000)にする関数
			void ClearUpperCounter(void);			//!< 上位カウンタをクリアする関数
		
		private:
			PCI6205C(const PCI6205C&) = delete;					//!< コピーコンストラクタ使用禁止
			const PCI6205C& operator=(const PCI6205C&) = delete;//!< 代入演算子使用禁止
			const unsigned int NUM_CH;				//!< 使用チャネル数
			const unsigned int ADDR_BASE0;			//!< 先頭アドレス
			const unsigned int ADDR_BASE1;			//!< 先頭アドレス
			const unsigned int ADDR_BASE2;			//!< 先頭アドレス
			const unsigned int ADDR_BASE3;			//!< 先頭アドレス
			unsigned int ADDR_COUNTER_LO[MAX_CH];	//!< 下位カウンタアドレス
			unsigned int ADDR_COUNTER_MD[MAX_CH];	//!< 中位カウンタアドレス
			unsigned int ADDR_COUNTER_HI[MAX_CH];	//!< 上位カウンタアドレス
			unsigned int ADDR_WRITEREG[MAX_CH];		//!< ライトレジスタアドレス
			unsigned int ADDR_STATUS[MAX_CH];		//!< ステータスレジスタアドレス
			uint8_t StatusReg[MAX_CH];				//!< ステータスレジスタ
			long CountPrev[MAX_CH];					//!< 前回のエンコーダデータ
			long UpperCount[MAX_CH];				//!< 上位カウンタ
			bool IsUsedZpulse;						//!< Zパルスを使用するか？
			
			void Settings(bool Mult4Mode);					// エンコーダカウンタの設定を行う関数
			long GetEncData(unsigned int ch);				// エンコーダカウント値を取得する関数
			void SetEncData(long count, unsigned int ch);	// エンコーダカウント値を設定する関数
			void SetCompData(long count, unsigned int ch);	// 比較レジスタにカウント値を設定する関数
			void ReadStatusReg(void);						// ステータスレジスタを読み込む関数  注意！：この関数は制御周期中に１回だけ呼び出すこと。
			bool GetDirectionFlag(unsigned int ch);			// 回転方向を返す関数  true=正方向, false=負方向
			bool GetCarryBorrowFlag(unsigned int ch);		// キャリー/ボローフラグを読み込んで返す関数  注意！：この関数を呼ぶ前に ReadStatusReg() を呼び出さないと最新の状態は得られない。
			bool GetEquDetectionFlag(unsigned int ch);		// カウント値比較検出フラグを読み込んで返す関数
			void LatchCounter(void);						// エンコーダカウンタを読み出しレジスタに移動(ラッチ)
			void SelectCounter(unsigned int ch);			// カウンタを選択する関数
			void SelectComparator(unsigned int ch);			// 比較レジスタを選択する関数
			static long ProcessCounterSign(long count);		// カウンタの符号処理
			long ProcessUpperCount(long count, unsigned int ch); // 上位カウンタの処理
			static unsigned long IIIbyteCat(unsigned short High, unsigned short Middle, unsigned short Low);// 上位、中位、下位に分かれている各々1バイトのデータを、3バイトのデータに結合する関数
			static uint8_t Get4byteHiHi(uint32_t in);	// 4byteデータの上位1byteを抽出して出力
			static uint8_t Get4byteHiLo(uint32_t in);	// 4byteデータの中上位1byteを抽出して出力
			static uint8_t Get4byteLoHi(uint32_t in);	// 4byteデータの中下位1byteを抽出して出力
			static uint8_t Get4byteLoLo(uint32_t in);	// 4byteデータの下位1byteを抽出して出力
			static bool GetBitStatus(uint8_t in, uint8_t bit);	// ビットの状態を返す関数  in : 変数, bit : ビット位置, 戻り値 true=1, false=0
	};
}

#endif

