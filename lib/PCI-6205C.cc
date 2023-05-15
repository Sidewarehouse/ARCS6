//! @file PCI-6205C.cc
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

#include "PCI-6205C.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

//! @brief コンストラクタ(ENC初期化＆設定)
//! @param[in] Addr0 ベースアドレス0
//! @param[in] Addr1 ベースアドレス1
//! @param[in] Addr2 ベースアドレス2
//! @param[in] Addr3 ベースアドレス3
//! @param[in] NumOfChannel 使用するチャネル数
//! @param[in] Mult4Mode 4逓倍モードを使用するか否か true = 使用，false = 非使用
PCI6205C::PCI6205C(const unsigned int Addr0, const unsigned int Addr1, const unsigned int Addr2, const unsigned int Addr3,
					unsigned int NumOfChannel, bool Mult4Mode)
	: NUM_CH(NumOfChannel), ADDR_BASE0(Addr0), ADDR_BASE1(Addr1), ADDR_BASE2(Addr2), ADDR_BASE3(Addr3),
	  ADDR_COUNTER_LO{}, ADDR_COUNTER_MD{}, ADDR_COUNTER_HI{}, ADDR_WRITEREG{}, ADDR_STATUS{},
	  StatusReg(), CountPrev{}, UpperCount{}, IsUsedZpulse(false)
{
	PassedLog();
	iopl(3);	// I/O全アドレス空間にアクセス許可
	
	// カウンタアドレスの計算
	ADDR_COUNTER_LO[0] = ADDR_BASE0 + 0x00;	// CH1下位カウンタアドレス
	ADDR_COUNTER_MD[0] = ADDR_BASE0 + 0x01;	// CH1中位カウンタアドレス
	ADDR_COUNTER_HI[0] = ADDR_BASE0 + 0x02;	// CH1上位カウンタアドレス
	ADDR_COUNTER_LO[1] = ADDR_BASE0 + 0x10;	// CH2下位カウンタアドレス
	ADDR_COUNTER_MD[1] = ADDR_BASE0 + 0x11;	// CH2中位カウンタアドレス
	ADDR_COUNTER_HI[1] = ADDR_BASE0 + 0x12;	// CH2上位カウンタアドレス
	ADDR_COUNTER_LO[2] = ADDR_BASE1 + 0x00;	// CH3下位カウンタアドレス
	ADDR_COUNTER_MD[2] = ADDR_BASE1 + 0x01;	// CH3中位カウンタアドレス
	ADDR_COUNTER_HI[2] = ADDR_BASE1 + 0x02;	// CH3上位カウンタアドレス
	ADDR_COUNTER_LO[3] = ADDR_BASE1 + 0x10;	// CH4下位カウンタアドレス
	ADDR_COUNTER_MD[3] = ADDR_BASE1 + 0x11;	// CH4中位カウンタアドレス
	ADDR_COUNTER_HI[3] = ADDR_BASE1 + 0x12;	// CH4上位カウンタアドレス
	ADDR_COUNTER_LO[4] = ADDR_BASE2 + 0x00;	// CH5下位カウンタアドレス
	ADDR_COUNTER_MD[4] = ADDR_BASE2 + 0x01;	// CH5中位カウンタアドレス
	ADDR_COUNTER_HI[4] = ADDR_BASE2 + 0x02;	// CH5上位カウンタアドレス
	ADDR_COUNTER_LO[5] = ADDR_BASE2 + 0x10;	// CH6下位カウンタアドレス
	ADDR_COUNTER_MD[5] = ADDR_BASE2 + 0x11;	// CH6中位カウンタアドレス
	ADDR_COUNTER_HI[5] = ADDR_BASE2 + 0x12;	// CH6上位カウンタアドレス
	ADDR_COUNTER_LO[6] = ADDR_BASE3 + 0x00;	// CH7下位カウンタアドレス
	ADDR_COUNTER_MD[6] = ADDR_BASE3 + 0x01;	// CH7中位カウンタアドレス
	ADDR_COUNTER_HI[6] = ADDR_BASE3 + 0x02;	// CH7上位カウンタアドレス
	ADDR_COUNTER_LO[7] = ADDR_BASE3 + 0x10;	// CH8下位カウンタアドレス
	ADDR_COUNTER_MD[7] = ADDR_BASE3 + 0x11;	// CH8中位カウンタアドレス
	ADDR_COUNTER_HI[7] = ADDR_BASE3 + 0x12;	// CH8上位カウンタアドレス
	
	// ライトレジスタアドレスの計算
	ADDR_WRITEREG[0] = ADDR_BASE0 + 0x05;	// CH1ライトレジスタアドレス
	ADDR_WRITEREG[1] = ADDR_BASE0 + 0x15;	// CH2ライトレジスタアドレス
	ADDR_WRITEREG[2] = ADDR_BASE1 + 0x05;	// CH3ライトレジスタアドレス
	ADDR_WRITEREG[3] = ADDR_BASE1 + 0x15;	// CH4ライトレジスタアドレス
	ADDR_WRITEREG[4] = ADDR_BASE2 + 0x05;	// CH5ライトレジスタアドレス
	ADDR_WRITEREG[5] = ADDR_BASE2 + 0x15;	// CH6ライトレジスタアドレス
	ADDR_WRITEREG[6] = ADDR_BASE3 + 0x05;	// CH7ライトレジスタアドレス
	ADDR_WRITEREG[7] = ADDR_BASE3 + 0x15;	// CH8ライトレジスタアドレス
	
	// ステータスレジスタアドレスの計算
	ADDR_STATUS[0] = ADDR_BASE0 + 0x06;	// CH1ステータスレジスタアドレス
	ADDR_STATUS[1] = ADDR_BASE0 + 0x16;	// CH2ステータスレジスタアドレス
	ADDR_STATUS[2] = ADDR_BASE1 + 0x06;	// CH3ステータスレジスタアドレス
	ADDR_STATUS[3] = ADDR_BASE1 + 0x16;	// CH4ステータスレジスタアドレス
	ADDR_STATUS[4] = ADDR_BASE2 + 0x06;	// CH5ステータスレジスタアドレス
	ADDR_STATUS[5] = ADDR_BASE2 + 0x16;	// CH6ステータスレジスタアドレス
	ADDR_STATUS[6] = ADDR_BASE3 + 0x06;	// CH7ステータスレジスタアドレス
	ADDR_STATUS[7] = ADDR_BASE3 + 0x16;	// CH8ステータスレジスタアドレス
	
	Settings(Mult4Mode);	// エンコーダ設定
	ClearAllCounter();		// カウンタクリア
}

//! @brief 空コンストラクタ
PCI6205C::PCI6205C()
	: NUM_CH(0), ADDR_BASE0(0), ADDR_BASE1(0), ADDR_BASE2(0), ADDR_BASE3(0),
	  ADDR_COUNTER_LO{}, ADDR_COUNTER_MD{}, ADDR_COUNTER_HI{}, ADDR_WRITEREG{}, ADDR_STATUS{},
	  StatusReg(), CountPrev{}, UpperCount{}, IsUsedZpulse(false)
{
	// 空即是色
}

//! @brief デストラクタ(ENC終了処理)
PCI6205C::~PCI6205C(){
	PassedLog();
	ClearAllCounter();
}

//! @brief エンコーダカウンタの設定を行う関数
//! @param[in] Mult4Mode 4逓倍モードを使用するか否か true = 使用，false = 非使用
void PCI6205C::Settings(bool Mult4Mode){
	if(Mult4Mode == false){
		outb(0x0C,ADDR_BASE0+0x04);	// CH1 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0C,ADDR_BASE0+0x14);	// CH2 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0C,ADDR_BASE1+0x04);	// CH3 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0C,ADDR_BASE1+0x14);	// CH4 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0C,ADDR_BASE2+0x04);	// CH5 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0C,ADDR_BASE2+0x14);	// CH6 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0C,ADDR_BASE3+0x04);	// CH7 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0C,ADDR_BASE3+0x14);	// CH8 位相差パルス 1逓倍 同期クリア 一致検出無効
		usleep(100);
	}else if(Mult4Mode == true){
		outb(0x0E,ADDR_BASE0+0x04);	// CH1 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0E,ADDR_BASE0+0x14);	// CH2 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0E,ADDR_BASE1+0x04);	// CH3 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0E,ADDR_BASE1+0x14);	// CH4 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0E,ADDR_BASE2+0x04);	// CH5 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0E,ADDR_BASE2+0x14);	// CH6 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0E,ADDR_BASE3+0x04);	// CH7 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
		outb(0x0E,ADDR_BASE3+0x14);	// CH8 位相差パルス 4逓倍 同期クリア 一致検出無効
		usleep(100);
	}
	ZpulseClear(true);	// Z相でカウンタクリア
	for(unsigned int i=0;i<NUM_CH;i++)SetCompData(0x000000, i);	// 比較レジスタを零に設定
}

//! @brief Z相クリア設定
//! @param[in] flag Z相パルスが来たらクリアするか否か true = クリアする，false = クリアしない
void PCI6205C::ZpulseClear(bool flag){
	IsUsedZpulse = flag;	// Z相パルス使用フラグを保持
	if(flag == true){
		outb(0x01,ADDR_BASE0+0x07);	// CH1 Z相でカウンタクリア
		outb(0x01,ADDR_BASE0+0x17);	// CH2 Z相でカウンタクリア
		outb(0x01,ADDR_BASE1+0x07);	// CH3 Z相でカウンタクリア
		outb(0x01,ADDR_BASE1+0x17);	// CH4 Z相でカウンタクリア
		outb(0x01,ADDR_BASE2+0x07);	// CH5 Z相でカウンタクリア
		outb(0x01,ADDR_BASE2+0x17);	// CH6 Z相でカウンタクリア
		outb(0x01,ADDR_BASE3+0x07);	// CH7 Z相でカウンタクリア
		outb(0x01,ADDR_BASE3+0x17);	// CH8 Z相でカウンタクリア
	}else{
		outb(0x00,ADDR_BASE0+0x07);	// CH1 カウンタクリアを無効
		outb(0x00,ADDR_BASE0+0x17);	// CH2 カウンタクリアを無効
		outb(0x00,ADDR_BASE1+0x07);	// CH3 カウンタクリアを無効
		outb(0x00,ADDR_BASE1+0x17);	// CH4 カウンタクリアを無効
		outb(0x00,ADDR_BASE2+0x07);	// CH5 カウンタクリアを無効
		outb(0x00,ADDR_BASE2+0x17);	// CH6 カウンタクリアを無効
		outb(0x00,ADDR_BASE3+0x07);	// CH7 カウンタクリアを無効
		outb(0x00,ADDR_BASE3+0x17);	// CH8 カウンタクリアを無効
	}
}

//! @brief エンコーダカウンタからカウント値を読み込む関数(生配列版)
//! @param[out] count エンコーダカウント値
void PCI6205C::GetCount(long count[MAX_CH]){
	unsigned int i;
	long buff[MAX_CH];	// カウント値読み込みバッファ
	LatchCounter();									// エンコーダカウンタをラッチしてから，
	for(i=0;i<NUM_CH;i++)buff[i] = GetEncData(i);	// エンコーダカウント値を取得
	if(IsUsedZpulse == true){
		for(i=0;i<NUM_CH;i++)buff[i] = ProcessCounterSign(buff[i]);	// Z相クリアが有効なときはカウンタの符号処理のみ実行
	}else{
		for(i=0;i<NUM_CH;i++)buff[i] = ProcessUpperCount(buff[i],i);// 逆に無効なときは上位カウンタの処理を実行 (Z相クリア使用して且つ上位カウントすると値がおかしくなる)
	}
	std::copy(buff, buff + MAX_CH, count);			// カウント値を返す
}

//! @brief エンコーダカウンタからカウント値を読み込む関数(std::array版)
//! @param[out] count エンコーダカウント値
void PCI6205C::GetCount(std::array<long, MAX_CH>& count){
	GetCount(count.data());
}

//! @brief エンコーダカウンタの値を零(0x000000)に設定，および上位カウンタをクリアする関数
void PCI6205C::ClearAllCounter(void){
	for(unsigned int i=0;i<NUM_CH;i++)SetEncData(0x000000, i);	// カウンタクリア
	ClearUpperCounter();	// 上位カウンタをクリア
}

//! @brief 上位カウンタをクリアする関数
void PCI6205C::ClearUpperCounter(void){
	std::fill(CountPrev, CountPrev + MAX_CH, 0);	// 前回のカウント値クリア
	std::fill(UpperCount, UpperCount + MAX_CH, 0);	// 上位カウンタクリア
}

//! @brief エンコーダカウント値を取得する関数
//! @param[in] ch チャネル番号(0から始まる)
long PCI6205C::GetEncData(unsigned int ch){
	uint8_t Hbuff, Mbuff, Lbuff;	// カウンタ値取得バッファ
	SelectCounter(ch);				// カウンタ選択
	Lbuff=inb(ADDR_COUNTER_LO[ch]);	// 下位 8bit データ取得
	Mbuff=inb(ADDR_COUNTER_MD[ch]);	// 中位 8bit データ取得
	Hbuff=inb(ADDR_COUNTER_HI[ch]);	// 上位 8bit データ取得
	return (long)IIIbyteCat(Hbuff,Mbuff,Lbuff);	// 上位、中位、下位データを24itのデータに結合する
}

//! @brief エンコーダカウント値を設定する関数
//! @param[in] count カウント値
//! @param[in] ch チャネル番号(0から始まる)
void PCI6205C::SetEncData(long count, unsigned int ch){
	SelectCounter(ch);								// カウンタ選択
	outb(Get4byteLoLo(count), ADDR_COUNTER_LO[ch]);	// 下位カウンタへ書き込み
	outb(Get4byteLoHi(count), ADDR_COUNTER_MD[ch]);	// 中位カウンタへ書き込み
	outb(Get4byteHiLo(count), ADDR_COUNTER_HI[ch]);	// 上位カウンタへ書き込み
}

//! @brief 比較レジスタにカウント値を設定する関数
//! @param[in] ch チャネル番号(0から始まる)
void PCI6205C::SetCompData(long count, unsigned int ch){
	SelectComparator(ch);	// 比較レジスタ選択
	outb(Get4byteLoLo(count), ADDR_COUNTER_LO[ch]);	// 下位カウンタへ書き込み
	outb(Get4byteLoHi(count), ADDR_COUNTER_MD[ch]);	// 中位カウンタへ書き込み
	outb(Get4byteHiLo(count), ADDR_COUNTER_HI[ch]);	// 上位カウンタへ書き込み
}

//! @brief ステータスレジスタを読み込む関数
//!        注意！：この関数は制御周期中に１回だけ呼び出すこと。
void PCI6205C::ReadStatusReg(void){
	for(unsigned int i=0;i<NUM_CH;i++)StatusReg[i] = inb(ADDR_STATUS[i]);	// ステータスレジスタ破壊読み込み
}

//! @brief 回転方向を返す関数
//! @param[in] ch チャネル番号(0から始まる)
//! @return 回転方向  true=正方向, false=負方向
bool PCI6205C::GetDirectionFlag(unsigned int ch){
	return GetBitStatus(StatusReg[ch], 0x01);	// 回転方向フラグを返す
}

//! @brief キャリー/ボローフラグを読み込んで返す関数
//!        注意！：この関数を呼ぶ前に ReadStatusReg() を呼び出さないと最新の状態は得られない。
//! @param[in] ch チャネル番号(0から始まる)
//! @return キャリー/ボローフラグ
bool PCI6205C::GetCarryBorrowFlag(unsigned int ch){
	return GetBitStatus(StatusReg[ch], 0x02);	// キャリー/ボローフラグを返す
}

//! @brief カウント値比較検出フラグを読み込んで返す関数
//!        注意！：この関数を呼ぶ前に ReadStatusReg() を呼び出さないと最新の状態は得られない。
//! @param[in] ch チャネル番号(0から始まる)
//! @return カウンタ一致発生フラグ
bool PCI6205C::GetEquDetectionFlag(unsigned int ch){
	return GetBitStatus(StatusReg[ch], 0x10);	// カウンタ一致発生フラグを返す
}

//! @brief エンコーダカウンタを読み出しレジスタに移動(ラッチ)
//!        エンコーダデータラッチ(アドレスはステータスレジスタと同一)
void PCI6205C::LatchCounter(void){
	for(unsigned int i=0;i<NUM_CH;i++)outb(0x02,ADDR_STATUS[i]);	// カウンタ値を読み出しレジスタに移動
}

//! @brief カウンタを選択する関数
//! @param[in] ch チャネル番号(0から始まる)
void PCI6205C::SelectCounter(unsigned int ch){
	outb(0x00,ADDR_WRITEREG[ch]);	// ライトレジスタをカウンタ選択に設定
}

//! @brief 比較レジスタを選択する関数
//! @param[in] ch チャネル番号(0から始まる)
void PCI6205C::SelectComparator(unsigned int ch){
	outb(0x01,ADDR_WRITEREG[ch]);	// ライトレジスタを比較レジスタ選択に設定
}
//! @brief カウンタの符号処理
//! @count[in] count エンコーダカウンタの生データ
//! @return 符号が考慮されたカウント値
long PCI6205C::ProcessCounterSign(long count){
	long ret = count;
	
	if(0x800000 <= (unsigned long)ret){
		// カウンタの最上位ビットが立っていたら，
		ret = (long)((unsigned long)ret | (unsigned long)0xFFFFFFFFFF000000);	// それより上位はすべて1で埋める
	}
	
	return ret;
}

//! @brief 上位カウンタの処理
//! @param[in] count エンコーダ生データ
//! @param[in] ch チャネル番号(0から始まる)
//! @return 上位カウンタと下位カウンタが結合された値
long PCI6205C::ProcessUpperCount(long count, unsigned int ch){
	if( 0x0000000000FFF000 < CountPrev[ch] && count < 0x0000000000000FFF ){	// 24bitオーバーフローを検出したら
		UpperCount[ch]++;					// 上位カウンタをカウントアップ
	}
	if( CountPrev[ch] < 0x0000000000000FFF && 0x0000000000FFF000 < count ){	// 24bitアンダーフローを検出したら
		UpperCount[ch]--;					// 上位カウンタをカウントダウン
	}
	CountPrev[ch] = count;					// 次回のために今回のカウント値を記録
	return (UpperCount[ch] << 24) | count;	// 上位カウンタと下位カウンタを結合して返す
}

//! @brief 上位、中位、下位に分かれている各々1バイトのデータを、3バイトのデータに結合する
//!        例： 0xAB, 0xCD, 0xEF → 0xABCDEF
//! @param[in] High 上位データ
//! @param[in] Middle 中位データ
//! @param[in] Low 下位データ
//! @return 上位・中位・下位が結合されたデータ
unsigned long PCI6205C::IIIbyteCat(unsigned short High, unsigned short Middle, unsigned short Low){
	unsigned long Hbuff=0, Mbuff=0, Lbuff=0, buff=0;
	Hbuff=(unsigned long)High   & 0x000000FF;
	Mbuff=(unsigned long)Middle & 0x000000FF;
	Lbuff=(unsigned long)Low    & 0x000000FF;
	
	// Hbuffを8bit左シフトしてMbuffとORを取り、その結果をさらに8bit左シフトしLbuffとORを取る
	buff = Hbuff << 8;
	buff = buff | Mbuff;
	buff = buff << 8;
	buff = buff | Lbuff;
	return buff;
}

//! @brief 4byteデータの上位1byteを抽出して出力
//! @param[in] in 4バイトデータ
//! @return 上位1バイトデータ
uint8_t PCI6205C::Get4byteHiHi(uint32_t in){
	return (uint8_t)((uint32_t)0x000000FF & (uint32_t)(in >> 24));
}

//! @brief 4byteデータの中上位1byteを抽出して出力
//! @param[in] in 4バイトデータ
//! @return 中上位1バイトデータ
uint8_t PCI6205C::Get4byteHiLo(uint32_t in){
	return (uint8_t)((uint32_t)0x000000FF & (uint32_t)(in >> 16));
}

//! @brief 4byteデータの中下位1byteを抽出して出力
//! @param[in] in 4バイトデータ
//! @return 中下位1バイトデータ
uint8_t PCI6205C::Get4byteLoHi(uint32_t in){
	return (uint8_t)((uint32_t)0x000000FF & (uint32_t)(in >> 8));
}

//! @brief 4byteデータの下位1byteを抽出して出力
//! @param[in] in 4バイトデータ
//! @return 下位1バイトデータ
uint8_t PCI6205C::Get4byteLoLo(uint32_t in){
	return (uint8_t)((uint32_t)0x000000FF & (uint32_t)in);
}

//! @brief ビットの状態を返す関数
//! @param[in] in 変数
//! @param[in] bit ビット位置
//! @return ビットの状態，戻り値 true = 1, false = 0
bool PCI6205C::GetBitStatus(uint8_t in, uint8_t bit){
	bool ret;
	if((in & bit) == bit){	// マスクしてから比較
		ret = true;			// ビットが1
	}else{
		ret = false;		// ビットが0
	}
	return ret;
}

