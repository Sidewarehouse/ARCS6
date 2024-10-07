//! @file PCI-3340.cc
//! @brief PCI-3340入出力クラス
//!
//! Interface社製PCI-3340のための入出力機能を提供します。(仕様変更済み)
//!
//! @date 2024/10/07
//! @author Yokokura, Yuki & Akira YAMAGUCHI
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// x86_64系の場合のみ対応-ここから
#ifdef __x86_64__

#include <sys/io.h>
#include <unistd.h>
#include <algorithm>
#include "PCI-3340.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

//! @brief コンストラクタ(DAC初期化＆設定)
//! @param[in] Base0 ベースアドレス0
PCI3340::PCI3340(unsigned int Base0)
	: BaseAddr0(Base0)	// Baseアドレスの格納
{
	PassedLog();
	Settings();
	Zero();
}

//! @brief 空コンストラクタ
PCI3340::PCI3340()
	: BaseAddr0(0)	// Baseアドレスの格納
{
	// 色即是空
}

//! @brief デストラクタ(DAC終了処理)
PCI3340::~PCI3340(){
	PassedLog();
	Zero();
}

//! @brief 指定した電圧をDACから出力
//! @param[in] Vout 出力電圧 [V]
void PCI3340::SetVoltage(double Vout[MAX_CH]){
	unsigned short DACdata[MAX_CH]={0x8000};
	for(unsigned int i=0;i<MAX_CH;i++)DACdata[i]=VoltToDacData(Vout[i]);	// 整数データに変換
	Output(DACdata);	// DA出力
}

//! @brief DAC出力電圧[V]からDACの実際の整数値に変換する
//! @param[in] Vdac 出力電圧 [V]
//! @return DACバイナリデータ
unsigned short PCI3340::VoltToDacData(double Vdac){
	// DACの仕様に従い定数を変更せよ
	// DACの設定：±10V 16bit (+10V=65525 0V=32767 -10V=0)
	return (unsigned short)( Vdac/(20.0/65536.0) + 32768.0 );
}

//! @brief DACの設定を行う関数
void PCI3340::Settings(void){
	iopl(3);	// I/O全アドレス空間にアクセス許可
	//------- DAC 設定
	// DAC チャネル 1
	outb(0x00,BaseAddr0+0x07);	// CH1 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	// DAC チャネル 2
	outb(0x01,BaseAddr0+0x07);	// CH2 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	// DAC チャネル 3
	outb(0x02,BaseAddr0+0x07);	// CH3 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	// DAC チャネル 4
	outb(0x03,BaseAddr0+0x07);	// CH4 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	// DAC チャネル 5
	outb(0x04,BaseAddr0+0x07);	// CH5 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	// DAC チャネル 6
	outb(0x05,BaseAddr0+0x07);	// CH6 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	// DAC チャネル 7
	outb(0x06,BaseAddr0+0x07);	// CH7 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	// DAC チャネル 8
	outb(0x07,BaseAddr0+0x07);	// CH8 設定
	outb(0x03,BaseAddr0+0x06);	// 出力電圧範囲を±10Vに設定
	usleep(100);				// 設定後に必要な100μsの待機
	
	outb(0x03,BaseAddr0+0x05);	// 全チャネル同時出力設定
	
	return;
}

//! @brief DACから指定した電圧を出力する関数
//! @param[in] DACdata DACバイナリデータ
void PCI3340::Output(unsigned short DACdata[MAX_CH]){
	//------- DAC 出力
	// DAC チャネル 1
	outb(0x00,BaseAddr0+0x02);					// CH1 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[0]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[0]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	// DAC チャネル 2
	outb(0x01,BaseAddr0+0x02);					// CH2 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[1]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[1]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	// DAC チャネル 3
	outb(0x02,BaseAddr0+0x02);					// CH3 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[2]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[2]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	// DAC チャネル 4
	outb(0x03,BaseAddr0+0x02);					// CH4 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[3]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[3]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	// DAC チャネル 5
	outb(0x04,BaseAddr0+0x02);					// CH5 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[4]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[4]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	// DAC チャネル 6
	outb(0x05,BaseAddr0+0x02);					// CH6 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[5]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[5]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	// DAC チャネル 7
	outb(0x06,BaseAddr0+0x02);					// CH7 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[6]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[6]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	// DAC チャネル 8
	outb(0x07,BaseAddr0+0x02);					// CH8 設定		[上位0xFF 下位0xFF]→ +10V
	outb(IIbyteLo(DACdata[7]),BaseAddr0+0x00);	// DAC出力 下位	[上位0x80 下位0x00]→   0V
	outb(IIbyteHi(DACdata[7]),BaseAddr0+0x01);	// DAC出力 上位	[上位0x00 下位0x00]→ -10V
	
	outb(0x01,BaseAddr0+0x05);	// 全チャネル同時出力実行
	
	return;
}

//! @brief DACの出力電圧を 0[V] にする関数
void PCI3340::Zero(void){
	unsigned short DACdata[MAX_CH]={0};
	for(unsigned int i=0;i<MAX_CH;i++)DACdata[i]=0x8000;
	Output(DACdata);
	return;
}

//! @brief 2byteデータの上位1byteを抽出して出力
//! @param[in] 2バイトデータ
//! @return 上位1バイトデータ
unsigned short PCI3340::IIbyteHi(unsigned short in){
	return 0x00FF & (in >> 8);
}

//! @brief 2byteデータの下位1byteを抽出して出力
//! @param[in] 2バイトデータ
//! @return 下位1バイトデータ
unsigned short PCI3340::IIbyteLo(unsigned short in){
	return 0x00FF & in;
}

//! @brief 汎用出力ピン2chからValを出力
//! @param[in] Val 出力バイナリデータ
void PCI3340::SetDigitalOut(unsigned short Val){
	outb(0x0003 & Val,BaseAddr0+0x1E);		//下位2bitをマスクしてValを出力
}

//! @brief 汎用入力ピン2chの値を入力
//! @return 入力バイナリデータ
unsigned short PCI3340::GetDigitalIn(void){
	return 0x0003 & inb(BaseAddr0+0x1E);	//下位2bitをマスクして入力
}

#endif
// x86_64系の場合のみ対応-ここまで
