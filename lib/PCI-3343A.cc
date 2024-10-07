//! @file PCI-3343A.cc
//! @brief PCI-3343A入出力クラス
//! Interface社製PCI-3343Aのための入出力機能を提供します。
//! 
//! @date 2024/10/07
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

// x86_64系の場合のみ対応-ここから
#ifdef __x86_64__

#include <sys/io.h>
#include <unistd.h>
#include <algorithm>
#include <stdint.h>
#include "PCI-3343A.hh"
#include "ARCSeventlog.hh"
#include "Limiter.hh"

using namespace ARCS;

//! @brief コンストラクタ(全チャネル使用する版)
//! @param [in] Addr ベースアドレス
PCI3343A::PCI3343A(unsigned int Addr)
  : ADDR_BASE(Addr),
	ADDR_DACDATA_LO(Addr + 0x00),
	ADDR_DACDATA_HI(Addr + 0x01),
	ADDR_CHSET(Addr + 0x02),
	ADDR_CONVMODE(Addr + 0x05),
	ADDR_OUTMODE(Addr + 0x1B),
	ADDR_DIO(Addr + 0x1E),
	ENA(0b00001111)
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
PCI3343A::PCI3343A(unsigned int Addr, uint8_t EnableCh)
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
PCI3343A::PCI3343A()
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
PCI3343A::~PCI3343A(){
	SetAllZero();			// 全チャネル零電圧出力設定
}

//! @brief 指定した電圧を出力する関数
//! @param [in] V1 CH1の電圧値[V]
//! @param [in] V2 CH2の電圧値[V]
//! @param [in] V3 CH3の電圧値[V]
//! @param [in] V4 CH4の電圧値[V]
void PCI3343A::SetVoltage(double V1, double V2, double V3, double V4){
	if( (ENA & 0b00000001) == 0b00000001 ){	// CH1がイネーブルのとき
		SelectCH(0);					// チャネル1選択
		SetDACdata(VoltToDacData(V1));	// DACデータ設定
	}
	if( (ENA & 0b00000010) == 0b00000010 ){	// CH2がイネーブルのとき
		SelectCH(1);					// チャネル2選択
		SetDACdata(VoltToDacData(V2));	// DACデータ設定
	}
	if( (ENA & 0b00000100) == 0b00000100 ){	// CH3がイネーブルのとき
		SelectCH(2);					// チャネル3選択
		SetDACdata(VoltToDacData(V3));	// DACデータ設定
	}
	if( (ENA & 0b00001000) == 0b00001000 ){	// CH4がイネーブルのとき
		SelectCH(3);					// チャネル4選択
		SetDACdata(VoltToDacData(V4));	// DACデータ設定
	}
	ExecOutput();	// 全チャネル同時電圧更新
}

//! @brief 全チャネル同時出力イネーブル
//! @param [in] flag true = 同時出力有効, false = 同時出力無効
void PCI3343A::SetAllEnable(bool flag){
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
void PCI3343A::SetOutEnable(bool flag){
	if(flag == true){
		outb(0x0F, ADDR_OUTMODE);	// 電圧出力有効
	}else{
		outb(0x00, ADDR_OUTMODE);	// ハイインピーダンスモード
	}
}

//! @brief 全チャネル同時出力実行(電圧更新)
void PCI3343A::ExecOutput(void){
	outb(0x01, ADDR_CONVMODE);
}

//! @brief チャネル選択
//! @param [in] ch チャネル番号
void PCI3343A::SelectCH(unsigned int ch){
	outb((unsigned char)ch, ADDR_CHSET);
}

//! @brief DACデータをセットする関数
//! @param [in] DACデータ
void PCI3343A::SetDACdata(uint16_t data){
	outb(Get2byteLo(data), ADDR_DACDATA_LO);	// 下位
	outb(Get2byteHi(data), ADDR_DACDATA_HI);	// 上位
}

//! @brief 全チャネル零電圧出力
void PCI3343A::SetAllZero(void){
	SelectCH(0);					// チャネル1選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(1);					// チャネル2選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(2);					// チャネル3選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	SelectCH(3);					// チャネル4選択
	SetDACdata(VoltToDacData(0));	// DACデータ設定
	ExecOutput();	// 全チャネル同時電圧更新
}

//! @brief DAC出力電圧[V]からDACの実際の整数値に変換する関数
//! @param [in] Vout 出力電圧
//! @return DACデータ
uint16_t PCI3343A::VoltToDacData(double Vout){
	// DACの設定：±10V 12bit (+10V=4096 0V=2048 -10V=0)
	return (uint16_t)( Limiter(Vout,10)*(4096.0/20.0) + 2048.0 );
}

//! @brief 2byteデータの上位1byteを抽出して出力
//! @param [in] 2byteデータ
//! @return 上位1byte
uint8_t PCI3343A::Get2byteHi(uint16_t in){
	return (uint8_t)((uint16_t)0x00FF & (uint16_t)(in >> 8));
}

//! @brief 2byteデータの下位1byteを抽出して出力
//! @param [in] 2byteデータ
//! @return 下位1byte
uint8_t PCI3343A::Get2byteLo(uint16_t in){
	return (uint8_t)((uint16_t)0x00FF & (uint16_t)in);
}

#endif
// x86_64系の場合のみ対応-ここまで
