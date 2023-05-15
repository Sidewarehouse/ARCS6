//! @file PCI-2826.cc
//! @brief PCI-2826入出力クラス
//! Interface社製PCI-2826のための入出力機能を提供します。
//! @date 2020/06/09
//! @author Yuto KOBAYASHI
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#include <sys/io.h>
#include <unistd.h>
#include <algorithm>
#include <stdint.h>
#include "PCI-2826CV.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;


//! @brief コンストラクタ
//! @param [in] Addr アドレス
PCI2826CV::PCI2826CV(unsigned int Addr)
  : InputData(0),
    OutputData(0),
    ADDR_BASE(Addr),
    ADDR_DIN0(Addr + 0x00),
    ADDR_DIN1(Addr + 0x01),
    ADDR_DIN2(Addr + 0x02),
    ADDR_DIN3(Addr + 0x03),
    ADDR_DOUT0(Addr + 0x00),
    ADDR_DOUT1(Addr + 0x01),
    ADDR_DOUT2(Addr + 0x02),
    ADDR_DOUT3(Addr + 0x03),
    GET_DATA_MASK(0xffffffff),
    SET_DATA_MASK(0xffffffff)
{
    PassedLog();
    iopl(3);        // I/O全アドレス空間にアクセス許可
    PCI2826CV::SetData(0);
}
//! @brief コンストラクタ
//! @param [in] Addr アドレス
//! @param [in] GetDataMask 入力用データマスク
//! @param [in] SetDataMask 出力用データマスク
PCI2826CV::PCI2826CV(unsigned int Addr, uint32_t GetDataMask, uint32_t SetDataMask)
  : InputData(0),
    OutputData(0),
    ADDR_BASE(Addr),
    ADDR_DIN0(Addr + 0x00),
    ADDR_DIN1(Addr + 0x01),
    ADDR_DIN2(Addr + 0x02),
    ADDR_DIN3(Addr + 0x03),
    ADDR_DOUT0(Addr + 0x00),
    ADDR_DOUT1(Addr + 0x01),
    ADDR_DOUT2(Addr + 0x02),
    ADDR_DOUT3(Addr + 0x03),
    GET_DATA_MASK(GetDataMask),
    SET_DATA_MASK(SetDataMask)
{
    PassedLog();
    iopl(3);        // I/O全アドレス空間にアクセス許可
    PCI2826CV::SetData(0);
}
//! @brief 空コンストラクタ
PCI2826CV::PCI2826CV()
 :  InputData(0),
    OutputData(0),
    ADDR_BASE(0),
    ADDR_DIN0(0),
    ADDR_DIN1(0),
    ADDR_DIN2(0),
    ADDR_DIN3(0),
    ADDR_DOUT0(0),
    ADDR_DOUT1(0),
    ADDR_DOUT2(0),
    ADDR_DOUT3(0),
    GET_DATA_MASK(0),
    SET_DATA_MASK(0)
{
    // 虚無
}
//! @brief デストラクタ
PCI2826CV::~PCI2826CV(){
    PassedLog();
}

//---------------------------//
// IOポート入力読み取り関数
//---------------------------//
//! @brief　全入力の読み取り
//! @param [out] デジタル読み出し値
uint32_t PCI2826CV::GetData(void){
    unsigned int buff[4];
    unsigned int data;
    // 各ポートからデータを取得
    iopl(3);	// なぜか必要（これがないとoutbでセグフォ、CTRL_LOOP内では初回のみで良い）
    buff[0] = inb(ADDR_DIN0);
    buff[1] = inb(ADDR_DIN1);
    buff[2] = inb(ADDR_DIN2);
    buff[3] = inb(ADDR_DIN3);
    // データをまとめてひとつの変数にする
    data = (unsigned int)((buff[3]&0xff)<<24 | (buff[2]&0xff)<<16 | (buff[1]&0xff)<<8 | (buff[0]&0xff));
    // マスク処理
    data = data & GET_DATA_MASK;
    InputData = data;
    return data;
}

//! @brief　あるポート入力の読み取り
//! @param [in] PortNum ポート番号（0-3）
//! @param [out] デジタル読み出し値
uint8_t PCI2826CV::GetData(uint8_t PortNum){
    uint8_t buff;
    uint8_t data;
    // 各ポートからデータを取得
    iopl(3);	// なぜか必要（これがないとoutbでセグフォ、CTRL_LOOP内では初回のみで良い）
    buff = inb(ADDR_DIN0 + PortNum);
    // マスク処理
    data = buff & ((GET_DATA_MASK>>(PortNum*8))&0xff);
    InputData = DataUpdate32((uint32_t)data<<(PortNum*8), InputData, 0x000000ff<<(PortNum*8));
    return data;
}

//! @brief　あるポート・ビット入力の読み取り
//! @param [in] PortNum ポート番号  (0-3）
//! @param [in] BitNum ビット番号   (0-7)
//! @param [out] デジタル読み出し値
bool PCI2826CV::GetData(uint8_t PortNum, uint8_t BitNum){
    unsigned int buff;
    unsigned int data;
    // 各ポートからデータを取得
    iopl(3);	// なぜか必要（これがないとoutbでセグフォ、CTRL_LOOP内では初回のみで良い）
    buff = inb(ADDR_BASE + PortNum);
    // マスク処理
    data = (buff>>(BitNum)) & ((GET_DATA_MASK>>(PortNum*8+BitNum))&0x01);
    InputData = DataUpdate32((uint32_t)data<<(PortNum*8+BitNum), InputData, 0x01<<(PortNum*8+BitNum));
    return data;
}


//---------------------------//
// IOポート出力書き込み関数　（立ち上がりmin 5us, max 65us）
//---------------------------//
//! @brief 全出力の書き込み
//! @param [in] デジタル書き込み値
void PCI2826CV::SetData(uint32_t Data){
    uint8_t buff[4];        // 各ポートのデータ格納
    uint32_t MaskedData;    // マスキングされたData格納

    //データのマスク
    MaskedData = Data & SET_DATA_MASK;
    // データを分割して各ポート用データに格納
    buff[0] = (MaskedData>>0)&0xff;
    buff[1] = (MaskedData>>8)&0xff;
    buff[2] = (MaskedData>>16)&0xff;
    buff[3] = (MaskedData>>24)&0xff;
    
    // 各ポートに出力
    iopl(3);	// なぜか必要（これがないとoutbでセグフォ、CTRL_LOOP内では初回のみで良い）
    outb(buff[0], ADDR_DOUT0);
    outb(buff[1], ADDR_DOUT1);
    outb(buff[2], ADDR_DOUT2);
    outb(buff[3], ADDR_DOUT3);
    // 現在の出力データを保持
    OutputData = MaskedData;
    return;
}
//! @brief あるポートの出力の書き込み
//! @param [in] デジタル書き込み値
//! @param [in] PortNum ポート番号  (0-3）
void PCI2826CV::SetData(uint8_t Data, uint8_t PortNum){
    uint32_t buff;      //データ格納用

    //データのマスクし格納
    buff = Data & ((GET_DATA_MASK>>PortNum*8)&0xff);
    
    // ポートに出力
    iopl(3);	// なぜか必要（これがないとoutbでセグフォ、CTRL_LOOP内では初回のみで良い）
    outb(buff, ADDR_DOUT0 + PortNum);
    // 現在の出力データを保持
    OutputData = DataUpdate32(buff<<(PortNum*8), OutputData, 0x000000ff<<(PortNum*8));
    return;
}
//! @brief あるポート・ビットの出力の書き込み
//! @param [in] デジタル書き込み値
//! @param [in] PortNum ポート番号  (0-3）
//! @param [in] BitNum ビット番号   (0-7)
void PCI2826CV::SetData(uint8_t Data, uint8_t PortNum, uint8_t BitNum){
    uint8_t buff;      // データ格納用

    // データのマスクし格納
    buff = (Data!=0) & ((GET_DATA_MASK>>(PortNum*8+BitNum))&0x01);
    // BitNumだけ書き換え
    buff = DataUpdate8(buff<<BitNum, (OutputData>>(PortNum*8))&0xff, 0x01<<BitNum);
    
    // ポートに出力
    iopl(3);	// なぜか必要（これがないとoutbでセグフォ、CTRL_LOOP内では初回のみで良い）
    outb(buff, ADDR_DOUT0 + PortNum);
    // 現在の出力データを保持
    OutputData = DataUpdate32(((uint32_t)buff)<<(PortNum*8), OutputData, 0x000000ff<<(PortNum*8));
    return;
}

//! @brief 32ビットデータのあるビットを更新
//! @param [in] NewData     更新したいデータ
//! @param [in] OrgData     現在のデータ
//! @param [in] BitArray    更新したいビット配列（そのビットが1で更新、0で保持）
//! @param [out] 更新されたデータ
uint32_t PCI2826CV::DataUpdate32(uint32_t NewData, uint32_t OrgData, uint32_t BitArray){
    return ((NewData & BitArray) | (OrgData & ~BitArray));
}

//! @brief 8ビットデータのあるビットを更新
//! @param [in] NewData     更新したいデータ
//! @param [in] OrgData     現在のデータ
//! @param [in] BitArray    更新したいビット配列（そのビットが1で更新、0で保持）
//! @param [out] 更新されたデータ
uint8_t PCI2826CV::DataUpdate8(uint8_t NewData, uint8_t OrgData, uint8_t BitArray){
    return ((NewData & BitArray) | (OrgData & ~BitArray));
}