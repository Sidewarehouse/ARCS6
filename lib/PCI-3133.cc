//! @file PCI-3133.cc
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

#include <sys/io.h>
#include <unistd.h>
#include "PCI-3133.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

//! @brief コンストラクタ(ADC初期化＆設定)
PCI3133::PCI3133(const unsigned int Base)
	: BaseAddr(Base)	// ベースアドレスの格納
{
	PassedLog();
}

//! @brief 空コンストラクタ(ADC初期化＆設定)
PCI3133::PCI3133()
	: BaseAddr(0)	// ベースアドレスの格納
{
	// 空即是色
}

//! @brief デストラクタ(ADC終了処理)
PCI3133::~PCI3133(){
	PassedLog();
}

void PCI3133::Settings(void){
	// A/Dコンバータの設定を行う関数
	// ADCの設定を行う関数
	// あらかじめバイポーラ±10Vに設定を行っているが，
	// 変えたい場合は説明書に従ってボードのスイッチをいじること。
	
	iopl(3);					// I/O全アドレス空間にアクセス許可
	outb(0x80, BaseAddr+0x00);	// Ch1にセット(念のため)
	usleep(50);					// チャンネル切換え時間待機
	outb(0x01, BaseAddr+0x01);	// 差動モードに設定
	usleep(100);				// 念のための待機
	outb(0x00,BaseAddr+0x03);
	outb(0x40, BaseAddr+0x00);	// AD変換開始
	usleep(50);					// 念のための待機
}

unsigned short PCI3133::Input(void){
	// A/Dコンバータから電圧データを読み込む関数
	// AD変換を行う関数
	// 上位・下位8ビットを取得し，結合したものをデータとして出力する。
	// ただ，12ビットのADCなので上位は正確にいうと4ビットしかない。
	
	unsigned short  Hbuff=0x00, Lbuff=0x00;	// 変換値取得用バッファ
	unsigned short ADCdata;
	
	Lbuff=inb(BaseAddr+0x00);		// CH1 下位 8bit データ取得
	Hbuff=inb(BaseAddr+0x01);		// CH1 上位 8bit データ取得
	ADCdata=IIbyteCat(Hbuff,Lbuff);	// 上位、中位、下位データを24itのデータに結合する
	outb(0x40, BaseAddr+0x00);		// AD変換開始
	
	return ADCdata;
}

double PCI3133::AdcDataToVolt(unsigned short ADCdata){
	// 電圧データから電圧値[V]に変換する関数
	// AD変換が完了している場合，最上位ビットが1となっているので
	// 0x8000だけ余分に引いている。
	return (double)(ADCdata-0x8000)*10.0/4096.0-5.0;
}

unsigned long PCI3133::IIbyteCat(unsigned short High, unsigned short Low){
	// 上位、下位に分かれている各々1バイトのデータを、2バイトのデータに結合する
	// 例： 0xAB, 0xCD → 0xABCD
	unsigned short buff=0;
	// Highを8bit左シフトしてLowとORを取る
	buff = High << 8;
	buff = buff | Low;
	return buff;
}

double PCI3133::GetVoltage(void){
	// 電圧値[V]を取得
	return AdcDataToVolt(Input());
}

