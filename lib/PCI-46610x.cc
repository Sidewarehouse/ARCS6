//! @file PCI-46610x.cc
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

#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include "PCI-46610x.hh"
#include "ARCSassert.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

const std::string PCI46610x::STR_CR(1, 0x0d);	//!< キャリッジリターン文字のアスキーコードの定義
const std::string PCI46610x::STR_LF(1, 0x0a);	//!< ラインフィード文字のアスキーコードの定義

//! @brief コンストラクタ
//! @param[in] Addr ベースアドレス
PCI46610x::PCI46610x(unsigned long Addr)
	: fd(0), IOptr(),
	  ADDR_BASE(Addr)
{
	PassedLog();
	
	// I/Oレジスタをメモリにマップ
	fd = open("/dev/mem", O_RDWR | O_SYNC);
	arcs_assert(fd != -1);
	IOptr = (uint8_t*)mmap(NULL, IOREG_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, fd, ADDR_BASE);
	close(fd);
	arcs_assert(IOptr != (void*)(-1));
	
	// 送受信FIFOメモリと設定レジスタを初期化
	ClearFIFOMemory();
	ClearAllSettings();
}

//! @brief ムーブコンストラクタ
PCI46610x::PCI46610x(PCI46610x&& right)
	: fd(0), IOptr(right.IOptr),
	  ADDR_BASE(right.ADDR_BASE)
{
	right.IOptr = nullptr;	// ムーブ元はヌルポにしとく
}

//! @brief デストラクタ
PCI46610x::~PCI46610x(){
	// I/Oレジスタをアンマップ
	munmap(IOptr, IOREG_SIZE);
	PassedLog();
}

//! @brief 通信設定の関数
void PCI46610x::SetConfig(enum BaudRate rate, enum Wire wir, enum Parity par, enum StopBit stpbit, enum WordLength wlen){
	uint8_t buff = 0;
	
	// ボーレートの設定
	switch(rate){
		// 2Mbpsのとき
		case RATE_2Mbps:
			IOptr[ADDR_BASECLOCK] = 0b00000000;	// 基準クロック周波数を32MHzに設定
			IOptr[ADDR_PRESCALER_DLL] = 0x01;	// 分周器を1に設定
			IOptr[ADDR_PRESCALER_DLM] = 0x00;	// 分周器を1に設定
			break;
		// 921.6kbpsのとき
		case RATE_921_6kbps:
			IOptr[ADDR_BASECLOCK] = 0b01000000;	// 基準クロック周波数を14.7456MHzに設定
			IOptr[ADDR_PRESCALER_DLL] = 0x01;	// 分周器を1に設定
			IOptr[ADDR_PRESCALER_DLM] = 0x00;	// 分周器を1に設定
			break;
		default:
			arcs_assert(false);	// ここには来ない
			break;
	}
	
	// 2線式/4線式の設定
	buff = IOptr[ADDR_DRIVE_CTRL];
	switch(wir){
		// 2線式のとき
		case WIRE_2:
			buff &= 0b11101111;	// 元から設定されてるビットを変化させないためにANDをとって指定ビットだけクリア
			break;
		// 4線式のとき
		case WIRE_4:
			buff |= 0b00010000;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		default:
			arcs_assert(false);	// ここには来ない
			break;
	}
	IOptr[ADDR_DRIVE_CTRL] = buff;
	
	// パリティの設定
	buff = IOptr[ADDR_LINE_CTRL];
	switch(par){
		// パリティ無効
		case PARITY_DISABLE:
			buff &= 0b11110111;	// 元から設定されてるビットを変化させないためにANDをとって指定ビットだけクリア
			break;
		// 奇数パリティ
		case PARITY_ODD:
			buff &= 0b11000111;	// 元から設定されてるビットを変化させないためにANDをとって指定ビットだけクリア
			buff |= 0b00001000;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		// 偶数パリティ
		case PARITY_EVEN:
			buff &= 0b11000111;	// 元から設定されてるビットを変化させないためにANDをとって指定ビットだけクリア
			buff |= 0b00011000;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		default:
			arcs_assert(false);	// ここには来ない
			break;
	}
	IOptr[ADDR_LINE_CTRL] = buff;
	
	// ストップビットの設定
	buff = IOptr[ADDR_LINE_CTRL];
	switch(stpbit){
		// ストップビット1
		case STOPBIT_1:
			buff &= 0b11111011;	// 元から設定されてるビットを変化させないためにANDをとって指定ビットだけクリア
			break;
		// ストップビット2
		case STOPBIT_2:
			buff |= 0b00000100;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		default:
			arcs_assert(false);	// ここには来ない
			break;
	}
	IOptr[ADDR_LINE_CTRL] = buff;
	
	// ワード長の設定
	buff = IOptr[ADDR_LINE_CTRL];
	buff &= 0b11111100;	// 元から設定されてるビットを変化させないためにANDをとって指定ビットだけクリア
	switch(wlen){
		// ワード長5bit
		case WORDLENG_5:
			buff |= 0b00000000;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		// ワード長6bit
		case WORDLENG_6:
			buff |= 0b00000001;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		// ワード長7bit
		case WORDLENG_7:
			buff |= 0b00000010;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		// ワード長8bit
		case WORDLENG_8:
			buff |= 0b00000011;	// 元から設定されてるビットを変化させないためにORをとって指定ビットだけセット
			break;
		default:
			arcs_assert(false);	// ここには来ない
			break;
	}
	IOptr[ADDR_LINE_CTRL] = buff;
	
	// T信号とC信号の設定
	buff = IOptr[ADDR_DRIVE_CTRL];
	buff |= 0b00000111;	// T信号とC信号を強制有効
	IOptr[ADDR_DRIVE_CTRL] = buff;
	
	EventLogVar(IOptr[ADDR_BASECLOCK]);
	EventLogVar(IOptr[ADDR_PRESCALER_DLL]);
	EventLogVar(IOptr[ADDR_PRESCALER_DLM]);
	EventLogVar(IOptr[ADDR_LINE_CTRL]);
	EventLogVar(IOptr[ADDR_DRIVE_CTRL]);
}

//! @brief 送信データ1byte分をセットする関数
//! @param[in]	TxData 1byte送信データ
void PCI46610x::SetTransData(uint8_t TxData){
	IOptr[ADDR_HOLDREG] = TxData;
}

//! @brief 受信データ1byte分を取得する関数
//! @return	RxData 1byte受信データ
uint8_t PCI46610x::GetRecvData(void) const{
	return IOptr[ADDR_HOLDREG];
}

//! @brief 全受信データを取得する関数
//! @param[out]	RxData 受信データの文字列
void PCI46610x::GetAllRecvData(std::string& RxData) const{
	// 受信FIFOメモリから1文字ずつ読み出す
	while(1){
		if(GetReceiveCount() == 0) break;				// 受信FIFOメモリの中身が無くなったらループを抜ける
		std::string buff(1, GetRecvData());				// FIFOから1文字取り出す
		if(buff == STR_CR || buff == STR_LF) buff = " ";// もし受信文字がCRかLFだったらスペースに変更
		RxData += buff;									// 受信文字をどんどん結合していく
	}
}

//! @brief 受信カウント数を返す関数
//! @return 受信カウント数(＝バイト数)
unsigned int PCI46610x::GetReceiveCount(void) const{
	uint16_t Lbuff = IOptr[ADDR_RXFIFO_CNT_LO];	// 下位8bit取得
	uint16_t Hbuff = IOptr[ADDR_RXFIFO_CNT_HI];	// 上位8bit取得
	return (Hbuff << 8) | Lbuff;	// 上位と下位を結合して返す
}

//! @brief 送受信FIFOメモリをクリアする関数
void PCI46610x::ClearFIFOMemory(void){
	uint8_t buff = IOptr[ADDR_FIFO_CTRL];	// 現在のレジスタ値を読み込み
	buff |= 0b00000110;								// TxFIFOとRxFIFOのリセットビットのみを立てる
	IOptr[ADDR_FIFO_CTRL] = buff;			// リセット実行
}

//! @brief すべての設定レジスタをクリアする関数
void PCI46610x::ClearAllSettings(void){
	IOptr[ADDR_RESET_REG] = 0b00000001;		// リセット実行
	// リセット完了まで待機する
	while(IOptr[ADDR_RESET_REG] == 0b00000001){
		asm("nop");	// 最適化防止のためのNOP
	}
}

