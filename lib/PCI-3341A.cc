//! @file PCI-3341A.cc
//! @brief PCI-3341A入出力クラス
//! Interface社製PCI-3341Aのための入出力機能を提供します。
//! PCI-3343Aのチャンネル数を拡張
//! @date 2024/04/09
//! @author Yokokura, Yuki Kosaka, Kohki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// x86_64系の場合のみ対応-ここから
#ifdef __x86_64__

#include <sys/io.h>
#include <unistd.h>
#include <algorithm>
#include <stdint.h>
#include <array>
#include "ArcsMatrix.hh"
#include "PCI-3341A.hh"
#include "ARCSeventlog.hh"
#include "Limiter.hh"


using namespace ARCS;

//! @brief コンストラクタ(全チャネル使用する版)
//! @param [in] Addr ベースアドレス
PCI3341A::PCI3341A(unsigned int Addr)
  : ADDR_BASE(Addr),
	ADDR_DACDATA_LO(Addr + 0x00),
	ADDR_DACDATA_HI(Addr + 0x01),
	ADDR_CHSET(Addr + 0x02),
	ADDR_CONVMODE(Addr + 0x05),
	ADDR_OUTMODE(Addr + 0x1B),
	ADDR_DIO(Addr + 0x1E),
	ENA(0b11111111)
{
	PassedLog();
	iopl(3);			// I/O全アドレス空間にアクセス許可
	SetAllEnable(true);	// 全チャネル同時出力を許可
	SetAllZero();		// 全チャネル零電圧出力設定
	SetOutEnable(true);	// 全チャネル出力有効
}

//! @brief コンストラクタ(指定チャネルのみ使用する版)
//! @param [in] Addr ベースアドレス
//! @param [in]	EnableCh	有効にするチャンネル(例：0b00000010ならCH2が有効)
PCI3341A::PCI3341A(unsigned int Addr, uint8_t EnableCh)
  : ADDR_BASE(Addr),
	ADDR_DACDATA_LO(Addr + 0x00),
	ADDR_DACDATA_HI(Addr + 0x01),
	ADDR_CHSET(Addr + 0x02),
	ADDR_CONVMODE(Addr + 0x05),
	ADDR_OUTMODE(Addr + 0x1B),
	ADDR_DIO(Addr + 0x1E),
	ENA(EnableCh)
{
	PassedLog();
	iopl(3);			// I/O全アドレス空間にアクセス許可
	SetAllEnable(true);	// 全チャネル同時出力を許可
	SetAllZero();		// 全チャネル零電圧出力設定
	SetOutEnable(true);	// 全チャネル出力有効
}

//! @brief 空コンストラクタ
PCI3341A::PCI3341A()
  : ADDR_BASE(0),
	ADDR_DACDATA_LO(0),
	ADDR_DACDATA_HI(0),
	ADDR_CHSET(0),
	ADDR_CONVMODE(0),
	ADDR_OUTMODE(0),
	ADDR_DIO(0),
	ENA(0)
{
	// 色即是空
}

//! @brief デストラクタ
PCI3341A::~PCI3341A(){
	SetAllZero();			// 全チャネル零電圧出力設定
}

//! @brief 指定した電圧を出力する関数
//! @param [in] Vout[n] CH[n+1]の電圧値[V]
void PCI3341A::SetVoltage(const ArcsMat<PCI3341A::MAX_CH, 1>& Vout){
	if( (ENA & 0b00000001) == 0b00000001 ){	// CH1がイネーブルのとき
		SelectCH(0);					// チャネル1選択
		SetDACdata(VoltToDacData(Vout[1]));	// DACデータ設定
	}
	if( (ENA & 0b00000010) == 0b00000010 ){	// CH2がイネーブルのとき
		SelectCH(1);					// チャネル2選択
		SetDACdata(VoltToDacData(Vout[2]));	// DACデータ設定
	}
	if( (ENA & 0b00000100) == 0b00000100 ){	// CH3がイネーブルのとき
		SelectCH(2);					// チャネル3選択
		SetDACdata(VoltToDacData(Vout[3]));	// DACデータ設定
	}
	if( (ENA & 0b00001000) == 0b00001000 ){	// CH4がイネーブルのとき
		SelectCH(3);					// チャネル4選択
		SetDACdata(VoltToDacData(Vout[4]));	// DACデータ設定
	}
    	if( (ENA & 0b00010000) == 0b00010000 ){	// CH5がイネーブルのとき
		SelectCH(4);					// チャネル5選択
		SetDACdata(VoltToDacData(Vout[5]));	// DACデータ設定
	}
    	if( (ENA & 0b00100000) == 0b00100000 ){	// CH6がイネーブルのとき
		SelectCH(5);					// チャネル6選択
		SetDACdata(VoltToDacData(Vout[6]));	// DACデータ設定
	}
    	if( (ENA & 0b0100000) == 0b0100000 ){	// CH7がイネーブルのとき
		SelectCH(6);					// チャネル7選択
		SetDACdata(VoltToDacData(Vout[7]));	// DACデータ設定
	}
    	if( (ENA & 0b10000000) == 0b10000000 ){	// CH8がイネーブルのとき
		SelectCH(7);					// チャネル8選択
		SetDACdata(VoltToDacData(Vout[8]));	// DACデータ設定
	}
	ExecOutput();	// 全チャネル同時電圧更新
}

//! @brief 全チャネル同時出力イネーブル
//! @param [in] flag true = 同時出力有効, false = 同時出力無効
void PCI3341A::SetAllEnable(bool flag){
	if(flag == true){
		// 同時出力許可
		outb(0x03, ADDR_CONVMODE);
	}else{
		// 同時出力禁止
		outb(0x00, ADDR_CONVMODE);
	}
}

//! @brief 全チャネル出力有効 
//! @param [in] flag = true で電圧出力，flag = false でハイインピーダンス
void PCI3341A::SetOutEnable(bool flag){
	if(flag == true){
		outb(0xFF, ADDR_OUTMODE);	// 電圧出力有効
	}else{
		outb(0x00, ADDR_OUTMODE);	// ハイインピーダンスモード
	}
}

//! @brief 全チャネル同時出力実行(電圧更新)
void PCI3341A::ExecOutput(void){
	outb(0x01, ADDR_CONVMODE);
}

//! @brief チャネル選択
//! @param [in] ch チャネル番号
void PCI3341A::SelectCH(unsigned int ch){
	outb((unsigned char)ch, ADDR_CHSET);
}

//! @brief DACデータをセットする関数
//! @param [in] DACデータ
void PCI3341A::SetDACdata(uint16_t data){
	outb(Get2byteLo(data), ADDR_DACDATA_LO);	// 下位
	outb(Get2byteHi(data), ADDR_DACDATA_HI);	// 上位
}

//! @brief 全チャネル零電圧出力
void PCI3341A::SetAllZero(void){
	SelectCH(0);					// チャネル1選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(1);					// チャネル2選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(2);					// チャネル3選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(3);					// チャネル4選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(4);					// チャネル4選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(5);					// チャネル4選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(6);					// チャネル4選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(7);					// チャネル4選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	ExecOutput();	// 全チャネル同時電圧更新
    
}

//! @brief DAC出力電圧[V]からDACの実際の整数値に変換する関数
//! @param [in] Vout 出力電圧
//! @return DACデータ
uint16_t PCI3341A::VoltToDacData(double Vout){
	// DACの設定：±10V 12bit (+10V=4096 0V=2048 -10V=0)
	return (uint16_t)( Limiter(Vout,10)*(4096.0/20.0) + 2048.0 );
}

//! @brief 2byteデータの上位1byteを抽出して出力
//! @param [in] 2byteデータ
//! @return 上位1byte
uint8_t PCI3341A::Get2byteHi(uint16_t in){
	return (uint8_t)((uint16_t)0x00FF & (uint16_t)(in >> 8));
}

//! @brief 2byteデータの下位1byteを抽出して出力
//! @param [in] 2byteデータ
//! @return 下位1byte
uint8_t PCI3341A::Get2byteLo(uint16_t in){
	return (uint8_t)((uint16_t)0x00FF & (uint16_t)in);
}


#endif
// x86_64系の場合のみ対応-ここまで
