//! @file PCI-3180.cc
//! @brief PCI-3180入出力クラス
//! Interface社製PCI-3180のための入出力機能を提供します。
//! @date 2024/10/07
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// 処理系でインクルードファイルを変える
#ifdef __x86_64__
	// x86_64系の場合
    #include <sys/io.h>
#endif
#ifdef __ARM_ARCH
	// ARM系の場合
    // sys/ioのインクルードなし
#endif

#include <unistd.h>
#include <algorithm>
#include <stdint.h>
#include "PCI-3180.hh"
#include "ARCSassert.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

//! @brief コンストラクタ(全チャネルを使用する版)
//! @param [in] Addr アドレス
//! @param [in] Range 電圧レンジ
PCI3180::PCI3180(unsigned int Addr, RangeMode Range)
  : ADDR_BASE(Addr),
	ADDR_ADCDATA_LO(Addr + 0x00),
	ADDR_ADCDATA_HI(Addr + 0x01),
	ADDR_CHSET(Addr + 0x02),
	ADDR_BUSY(Addr + 0x03),
	ADDR_CONVST(Addr + 0x04),
	ADDR_RANGE(Addr + 0x06),
	ADDR_RANGE_CH(Addr + 0x07),
	ADDR_DIO(Addr + 0x1E),
	A(0),
	B(0),
	ENA(0b00001111)
{
	PassedLog();
	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		iopl(3);			// I/O全アドレス空間にアクセス許可
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif
	
	// 入力電圧レンジの設定
	CalcVoltConv(Range);	// 電圧換算の傾きと切片の計算
	SelectRangeCH(0);		// チャネル1選択
	SetInputRange(Range);	// 電圧レンジ設定
	SelectRangeCH(1);		// チャネル2選択
	SetInputRange(Range);	// 電圧レンジ設定
	SelectRangeCH(2);		// チャネル3選択
	SetInputRange(Range);	// 電圧レンジ設定
	SelectRangeCH(3);		// チャネル4選択
	SetInputRange(Range);	// 電圧レンジ設定
}

//! @brief コンストラクタ(指定したチャネルのみを使用する版)
//! @param [in] Addr アドレス
//! @param [in] Range 電圧レンジ
//! @param [in] EnableCh イネーブルチャネル(有効にしたいチャネルの選択)
PCI3180::PCI3180(unsigned int Addr, RangeMode Range, uint8_t EnableCh)
  : ADDR_BASE(Addr),
	ADDR_ADCDATA_LO(Addr + 0x00),
	ADDR_ADCDATA_HI(Addr + 0x01),
	ADDR_CHSET(Addr + 0x02),
	ADDR_BUSY(Addr + 0x03),
	ADDR_CONVST(Addr + 0x04),
	ADDR_RANGE(Addr + 0x06),
	ADDR_RANGE_CH(Addr + 0x07),
	ADDR_DIO(Addr + 0x1E),
	A(0),
	B(0),
	ENA(EnableCh)
{
	PassedLog();

	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		iopl(3);			// I/O全アドレス空間にアクセス許可
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif
	
	// 入力電圧レンジの設定
	CalcVoltConv(Range);	// 電圧換算の傾きと切片の計算
	if( (ENA & 0b00000001) == 0b00000001 ){	// CH1がイネーブルのとき
		SelectRangeCH(0);		// チャネル1選択
		SetInputRange(Range);	// 電圧レンジ設定
	}
	if( (ENA & 0b00000010) == 0b00000010 ){	// CH2がイネーブルのとき
		SelectRangeCH(1);		// チャネル2選択
		SetInputRange(Range);	// 電圧レンジ設定
	}
	if( (ENA & 0b00000100) == 0b00000100 ){	// CH3がイネーブルのとき
		SelectRangeCH(2);		// チャネル3選択
		SetInputRange(Range);	// 電圧レンジ設定
	}
	if( (ENA & 0b00001000) == 0b00001000 ){	// CH4がイネーブルのとき
		SelectRangeCH(3);		// チャネル4選択
		SetInputRange(Range);	// 電圧レンジ設定
	}
}

//! @brief 空コンストラクタ
PCI3180::PCI3180()
  : ADDR_BASE(0),
	ADDR_ADCDATA_LO(0),
	ADDR_ADCDATA_HI(0),
	ADDR_CHSET(0),
	ADDR_BUSY(0),
	ADDR_CONVST(0),
	ADDR_RANGE(0),
	ADDR_RANGE_CH(0),
	ADDR_DIO(0),
	A(0),
	B(0),
	ENA(0)
{
	// 色即是空
}

//! @brief デストラクタ
PCI3180::~PCI3180(){
	PassedLog();
}

//! @brief AD変換開始
void PCI3180::ConvStart(void){
	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		outb(0x00, ADDR_CONVST);
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif
}

//! @brief 変換待機信号の取得
//! @return true = 変換中, false = 変換終了
bool PCI3180::GetBusy(void){
	bool ret = true;

	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		// 念のため 10000000 でAND取ってマスクしてから評価
		if((0x80 & inb(ADDR_BUSY)) == 0x80){
			ret = false;	// 変換終了
		}else{
			ret = true;		// 変換中
		}
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif

	return ret;
}

//! @brief AD変換が完了するまでビジーウェイト(ブロッキング動作関数)
void PCI3180::WaitBusy(void){
	while(GetBusy()){
		asm("nop");		// 最適化防止のための "NOP"
	}
}

//! @brief チャネル選択
//! @param [in] ch チャネル番号
void PCI3180::SelectCH(unsigned int ch){
	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		outb((unsigned char)ch, ADDR_CHSET);
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif
}

//! @brief ADCデータを取得する関数
//! @return ADC変換データ
uint16_t PCI3180::GetADCdata(void){
	uint8_t hi = 0, lo = 0;

	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		lo = inb(ADDR_ADCDATA_LO);	// 下位データの取得
		hi = inb(ADDR_ADCDATA_HI);	// 上位データの取得
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif

	return Combine2byte(hi, lo);	// 上位と下位を結合して返す
}

//! @brief 電圧を取得する関数
//! @param [out] V1 CH1の電圧値[V]
//! @param [out] V2 CH2の電圧値[V]
//! @param [out] V3 CH3の電圧値[V]
//! @param [out] V4 CH4の電圧値[V]
void PCI3180::GetVoltage(double& V1, double& V2, double& V3, double& V4){
	if( (ENA & 0b00000001) == 0b00000001 ){	// CH1がイネーブルのとき
		SelectCH(0);	// チャネル1選択
		V1 = AdcDataToVolt(GetADCdata());
	}else{
		V1 = 0;
	}
	if( (ENA & 0b00000010) == 0b00000010 ){	// CH2がイネーブルのとき
		SelectCH(1);	// チャネル2選択
		V2 = AdcDataToVolt(GetADCdata());
	}else{
		V2 = 0;
	}
	if( (ENA & 0b00000100) == 0b00000100 ){	// CH3がイネーブルのとき
		SelectCH(2);	// チャネル3選択
		V3 = AdcDataToVolt(GetADCdata());
	}else{
		V3 = 0;
	}
	if( (ENA & 0b00001000) == 0b00001000 ){	// CH4がイネーブルのとき
		SelectCH(3);	// チャネル4選択
		V4 = AdcDataToVolt(GetADCdata());
	}else{
		V4 = 0;
	}
}

//! @brief ADCデータから電圧値[V]に変換する関数
//! @param [in] ADCデータ
double PCI3180::AdcDataToVolt(uint16_t data) const {
	return A*(double)data + B;
}

//! @brief 入力レンジ設定チャンネル選択
//! @param [in] ch チャネル番号
void PCI3180::SelectRangeCH(unsigned int ch){
	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		outb((unsigned char)ch, ADDR_RANGE_CH);
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif
}

//! @brief 入力レンジの設定
//! @param [in] 入力電圧範囲
void PCI3180::SetInputRange(RangeMode Range){
	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		switch(Range){
			case RANGE_B_10V:
				outb(0x01, ADDR_RANGE);	// ±10Vに設定
				break;
			case RANGE_B_5V:
				outb(0x00, ADDR_RANGE);	// ±5Vに設定
				break;
			case RANGE_U_10V:
				outb(0x05, ADDR_RANGE);	// 0～10Vに設定
				break;
			default:
				outb(0x01, ADDR_RANGE);	// デフォルトは±10Vに設定
				break;
		}
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif
}

//! @brief 電圧換算の傾きと切片の計算
//! @param [in] 入力電圧範囲
void PCI3180::CalcVoltConv(RangeMode Range){
	switch(Range){
		case RANGE_B_10V:
			A = 2.0*10.0/4096.0;// 電圧換算の傾きの計算
			B = -10;			// 電圧換算の切片の計算
			break;
		case RANGE_B_5V:
			A = 2.0*5.0/4096.0;	// 電圧換算の傾きの計算
			B = -5;				// 電圧換算の切片の計算
			break;
		case RANGE_U_10V:
			A = 10.0/4096.0;	// 電圧換算の傾きの計算
			B = 0;				// 電圧換算の切片の計算
			break;
		default:
			A = 2.0*10.0/4096.0;// 電圧換算の傾きの計算
			B = -10;			// 電圧換算の切片の計算
			break;
	}
}

//! @brief 上位、下位に分かれている各々1バイトのデータを、2バイトのデータに結合する
//! 例： 0xAB, 0xCD → 0xABCD
//! @param [in] High 上位8bit
//! @param [in] Low  下位8bit
//! @return 結合された16bitデータ
uint16_t PCI3180::Combine2byte(uint8_t High, uint8_t Low){
	uint16_t buff;
	// Highを8bit左シフトしてLowとORを取る
	buff = High << 8;
	buff = buff | Low;
	return buff;
}

//! ディジタル出力ポートを設定する関数
//! @param [in] Data ディジタル出力バイナリデータ(2ch分)
void PCI3180::SetDigitalOutput(const uint8_t& Data){
	// 処理系で挙動を変える
	#ifdef __x86_64__
	    // x86_64系の場合
		outb((unsigned char)Data, ADDR_DIO);
	#endif
	#ifdef __ARM_ARCH
		// ARM系の場合
		arcs_assert(false);	// ARM系は未対応
	#endif
}
