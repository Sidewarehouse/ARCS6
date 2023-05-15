//! @file RPi2GPIO.cc
//! @brief Raspberry Pi 2 GPIOコントローラ
//!
//! ラズベリーパイ2の汎用入出力の制御を行うクラス
//!
//! @date 2019/08/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.
//
// ＜I/Oポート設定レジスタ → IOptr[0～3]にマッピング＞ (数字の意味はGPIO番号の1桁目を表す)
// GPIO09～GPIO 0 → 0x3F200003～0x3F200000 → IOptr[0] : --999888 77766655 54443332 22111000
// GPIO19～GPIO10 → 0x3F200007～0x3F200004 → IOptr[1] : --999888 77766655 54443332 22111000
// GPIO29～GPIO20 → 0x3F20000B～0x3F200008 → IOptr[2] : --999888 77766655 54443332 22111000
// GPIO39～GPIO30 → 0x3F20000F～0x3F20000C → IOptr[3] : --999888 77766655 54443332 22111000
// GPIO49～GPIO40 → 0x3F200013～0x3F200010 → IOptr[4] : --999888 77766655 54443332 22111000
// GPIO53～GPIO50 → 0x3F200017～0x3F200014 → IOptr[5] : -------- -------- ----3332 22111000
// ※各々のI/Oピンは下記のように3ビットのデータに従って設定される
//   000 = 入力，001 = 出力，[010,011,111,110,101,100] = Alternate Function[5,4,3,2,1,0]
//
// ＜未使用レジスタ → IOptr[6,9,12]にマッピング＞
// 0x3F20001B～0x3F200018 → IOptr[6]
// 0x3F200027～0x3F200024 → IOptr[9]
// 0x3F200033～0x3F200030 → IOptr[12]
//
// ＜I/Oポートレジスタ(セット用) → IOptr[7～8]にマッピング＞ (数字の意味はGPIO番号の1桁目を表す)
// GPIO31～GPIO 0 → 0x3F20001F～0x3F20001C → IOptr[7] : 10987654 32109876 54321098 76543210
// GPIO53～GPIO32 → 0x3F200023～0x3F200020 → IOptr[8] : -------- --321098 76543210 98765432
//
// ＜I/Oポートレジスタ(クリア用) → IOptr[10～11]にマッピング＞ (数字の意味はGPIO番号の1桁目を表す)
// GPIO31～GPIO 0 → 0x3F20002B～0x3F200028 → IOptr[10]: 10987654 32109876 54321098 76543210
// GPIO53～GPIO32 → 0x3F20002F～0x3F20002C → IOptr[11]: -------- --321098 76543210 98765432
//
// ＜I/Oポートレジスタ(読み込み用) → IOptr[13～14]にマッピング＞ (数字の意味はGPIO番号の1桁目を表す)
// GPIO31～GPIO 0 → 0x3F200037～0x3F200034 → IOptr[13]: 10987654 32109876 54321098 76543210
// GPIO53～GPIO32 → 0x3F20003B～0x3F200038 → IOptr[14]: -------- --321098 76543210 98765432

#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include "RPi2GPIO.hh"

using namespace ARCS;

//! @brief コンストラクタ
RPi2GPIO::RPi2GPIO()
	: fd(0), IOptr()
{
	// I/Oポートをメモリにマップ
	fd = open("/dev/mem", O_RDWR | O_SYNC);
	IOptr = (uint32_t*)mmap(NULL, IOREG_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, fd, IOREG_ADDRESS);
	close(fd);
	
	SetConfig_AllOutput();	// とりあえず最初はGPIO39～0を出力モードに設定
}

//! @brief デストラクタ
RPi2GPIO::~RPi2GPIO(){
	munmap(IOptr, IOREG_SIZE);	// I/Oポートをアンマップ
}

//! @brief ムーブコンストラクタ
RPi2GPIO::RPi2GPIO(RPi2GPIO&& r)
	: fd(r.fd), IOptr(r.IOptr)
{
	IOptr = nullptr;
}

void RPi2GPIO::SetConfig(int port, IOmode mode){
	// 入出力を設定する関数 port : GPIO番号，mode : RPi2GPIO::IN or RPi2GPIO::OUT
	// イジりたいビットのみをとりあえずゼロにする
	// (やってることは、GPIO番号の1桁目の数の分だけ左にビットシフトしてNOTしてから、元のビット列とANDを取る)
	IOptr[(unsigned int)(port/10)] &= ~( 0b00000000000000000000000000000111 << (unsigned int)((port % 10)*3) );
	// 出力モードならば指定部分に 001 を書き込み、入力モードなら既に 000 になってるはずなので何もしない
	if(mode == OUT){
		IOptr[(unsigned int)(port/10)] |= 0b00000000000000000000000000000001 << (unsigned int)((port % 10)*3);
	}
}

void RPi2GPIO::SetConfig_GPIO9to0(uint32_t bitdata){
	// GPIO9～0の設定レジスタを一括設定する関数
	IOptr[0] = bitdata;
}

void RPi2GPIO::SetConfig_GPIO19to10(uint32_t bitdata){
	// GPIO19～10の設定レジスタを一括設定する関数
	IOptr[1] = bitdata;
}

void RPi2GPIO::SetConfig_GPIO29to20(uint32_t bitdata){
	// GPIO29～20の設定レジスタを一括設定する関数
	IOptr[2] = bitdata;
}

void RPi2GPIO::SetConfig_GPIO39to30(uint32_t bitdata){
	// GPIO39～30の設定レジスタを一括設定する関数
	IOptr[3] = bitdata;
}

void RPi2GPIO::SetConfig_GPIO49to40(uint32_t bitdata){
	// GPIO49～40の設定レジスタを一括設定する関数
	IOptr[4] = bitdata;
}

void RPi2GPIO::SetConfig_GPIO53to50(uint32_t bitdata){
	// GPIO53～50の設定レジスタを一括設定する関数
	IOptr[5] = bitdata;
}

void RPi2GPIO::SetConfig_AllOutput(void){
	// GPIO39～0を出力モードに設定する関数
	//                     --999888777666555444333222111000
	SetConfig_GPIO9to0(  0b00001001001001001001001001001001);
	SetConfig_GPIO19to10(0b00001001001001001001001001001001);
	SetConfig_GPIO29to20(0b00001001001001001001001001001001);
	SetConfig_GPIO39to30(0b00001001001001001001001001001001);
	// GPIO53～40に書き込むとフリーズするので、放っておく
}

void RPi2GPIO::SetBitDataLo(uint32_t bitdata){
	// GPIO31～0 からビットデータを一括出力する関数
	IOptr[7]  = bitdata;	// ビットセット
	IOptr[10] = ~bitdata;	// ビットクリア
	msync(IOptr, IOREG_SIZE, MS_ASYNC);	// ファイルを同期
}

void RPi2GPIO::SetBitDataHi(uint32_t bitdata){
	// GPIO53～32 からビットデータを一括出力する関数
	IOptr[8]  = 0b00000000000111111111111111111111 & bitdata;	// 上位はマスクしてからビットセット
	IOptr[11] = 0b00000000000111111111111111111111 & ~bitdata;	// 上位はマスクしてからビットクリア
	msync(IOptr, IOREG_SIZE, MS_ASYNC);	// ファイルを同期
}

void RPi2GPIO::BitSet(unsigned int port){
	// 指定したGPIOを 1 にする関数
	if(port < 32){
		IOptr[7] = 0x00000001 << port;	// GPIO番号の分だけ左シフト
	}else if(32 <= port && port < 54){
		IOptr[8] = 0x00000001 << (port - 32);	// GPIO32以上のときはアドレスを1増加させて書き込む
	}
	msync(IOptr, IOREG_SIZE, MS_ASYNC);	// ファイルを同期
}

void RPi2GPIO::BitClear(unsigned int port){
	// 指定したGPIOを 0 にする関数
	if(port < 32){
		IOptr[10] = 0x00000001 << port;	// GPIO番号の分だけ左シフト
	}else if(32 <= port && port < 54){
		IOptr[11] = 0x00000001 << (port - 32);	// GPIO32以上のときはアドレスを1増加させて書き込む
	}
	msync(IOptr, IOREG_SIZE, MS_ASYNC);	// ファイルを同期
}

uint32_t RPi2GPIO::GetBitDataLo(void) const {
	// GPIO31～0 からデータを一括入力する関数
	msync(IOptr, IOREG_SIZE, MS_SYNC);	// ファイルを同期
	return (uint32_t)IOptr[13];
}

uint32_t RPi2GPIO::GetBitDataHi(void) const {
	// GPIO53～32 からデータを一括入力する関数
	msync(IOptr, IOREG_SIZE, MS_SYNC);	// ファイルを同期
	return (uint32_t)(0b00000000000111111111111111111111 & IOptr[14]);	// 上位はマスクしてから返す
}

bool RPi2GPIO::BitGet(unsigned int port) const {
	// 指定したGPIOからビットの状態を取得する関数
	msync(IOptr, IOREG_SIZE, MS_SYNC);	// ファイルを同期
	if(port < 32){
		return (bool)(0x00000001 & ((uint32_t)IOptr[13] >> port));	// GPIO番号の分だけ右シフトしてから上位をマスクして返す
	}else if(32 <= port && port < 54){
		return (bool)(0x00000001 & ((uint32_t)IOptr[14] >> (port - 32) ));	// GPIO32以上のときはアドレスを1増加させて読み込む
	}
	return false;
}

void RPi2GPIO::SetACTLED(bool onoff){
	// ACT LED の点灯を制御する関数 onoff=trueでオン、onff=falseでオフ
	if(onoff == true){
		BitSet(47);
	}else{
		BitClear(47);
	}
}

void RPi2GPIO::SetPWRLED(bool onoff){
	// PWR LED の点灯を制御する関数 onoff=trueでオン、onff=falseでオフ
	if(onoff == true){
		BitSet(35);
	}else{
		BitClear(35);
	}
}

