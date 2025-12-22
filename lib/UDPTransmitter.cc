//! @file UDPTransmitter.cc
//! @brief UDP送信器クラス
//!
//! 指定したIPアドレスおよびPort番号へUDPで配列変数値を飛ばすクラス
//! 後々、テンプレートに書き換え予定
//!
//! @date 2025/12/22
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2025 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include "UDPTransmitter.hh"

using namespace ARCS;

UDPTransmitter::UDPTransmitter(const char* IPaddress, const unsigned short PortNum, VarType TypeOfVar, unsigned int NumOfVar)
	// コンストラクタ (関数オーバーロード: 宛先のみ指定する場合)
	// IPaddress_dest : 宛先IPアドレス, PortNum_dest : 宛先Port番号
	// TypeOfVar : 送信する変数の型(ヘッダの定義を参照せよ), NumOfVar : 送信する変数の数
: sock(0), addr_dest(), addr_dept(), VariableType(TypeOfVar), TransData(0), TransBytes(0)
{
	sock = socket(AF_INET, SOCK_DGRAM, 0);	// ソケット生成
	
	// 宛先のIPとポートを指定
	memset(&addr_dest, 0, sizeof(addr_dest));			// 初期化
	addr_dest.sin_family = AF_INET;						// ファミリ指定
	addr_dest.sin_port = htons(PortNum);				// Port番号の設定
	addr_dest.sin_addr.s_addr = inet_addr(IPaddress);	// IPアドレスの設定
	
	TransBytes = CalcTransBytes(TypeOfVar, NumOfVar);	// 送信するデータの型からバイト数を計算
	TransData = new char[TransBytes];					// 送信バイナリデータ用の配列確保
}

UDPTransmitter::UDPTransmitter(const char* IPaddress_dept, const unsigned short PortNum_dept,
	const char* IPaddress_dest, const unsigned short PortNum_dest, VarType TypeOfVar, unsigned int NumOfVar)
	// コンストラクタ (関数オーバーロード: 発信元を指定する場合)
	// IPaddress_dept : 発信元IPアドレス, PortNum_dept : 発信元Port番号
	// IPaddress_dest : 宛先IPアドレス, PortNum_dest : 宛先Port番号
	// TypeOfVar : 送信する変数の型(ヘッダの定義を参照せよ), NumOfVar : 送信する変数の数
: sock(0), addr_dest(), addr_dept(), VariableType(TypeOfVar), TransData(0), TransBytes(0)
{
	sock = socket(AF_INET, SOCK_DGRAM, 0);	// ソケット生成
	
	// 発信元のIPとポートを指定
	memset(&addr_dept, 0, sizeof(addr_dept));				// 初期化
	addr_dept.sin_family = AF_INET;							// ファミリ指定
	addr_dept.sin_port = htons(PortNum_dept);				// Port番号の設定
	addr_dept.sin_addr.s_addr = inet_addr(IPaddress_dept);	// IPアドレスの設定
	
	// 宛先のIPとポートを指定
	memset(&addr_dest, 0, sizeof(addr_dest));				// 初期化
	addr_dest.sin_family = AF_INET;							// ファミリ指定
	addr_dest.sin_port = htons(PortNum_dest);				// Port番号の設定
	addr_dest.sin_addr.s_addr = inet_addr(IPaddress_dest);	// IPアドレスの設定
	
	// 送信元指定のためのバインド
	bind(sock, (struct sockaddr *)&addr_dept, sizeof(addr_dept));
	
	TransBytes = CalcTransBytes(TypeOfVar, NumOfVar);	// 送信するデータの型からバイト数を計算
	TransData = new char[TransBytes];					// 送信バイナリデータ用の配列確保
}

UDPTransmitter::~UDPTransmitter(){
	// デストラクタ
	// ソケット破棄
	close(sock);
	delete[] TransData;	// 送信バイナリデータ用配列の破棄
}

volatile bool UDPTransmitter::SendTo(void){
	// 実際に送信を実行する関数
	if(sendto(sock, TransData, TransBytes, 0, (struct sockaddr *)&addr_dest, sizeof(addr_dest))==-1){
		return false;	// 送信に失敗した場合は false を返す
	};
	return true;
}

size_t UDPTransmitter::CalcTransBytes(VarType TypeOfVar, unsigned int NumOfVar){
	// 送信するデータの型からバイト数を計算する関数 (そのうち気が向いたらテンプレートで書きなおすかも)
	// (わざわざこんなことをしているのは、このクラスを使う人がバイト数を意識しないで済むようにするため)
	
	size_t ret = 0;
	
	switch(TypeOfVar){
		case UDP_INT16_T:
			ret = sizeof(int16_t)*NumOfVar;	// 送信バイナリデータのバイト数
			break;
		case UDP_INT32_T:
			ret = sizeof(int32_t)*NumOfVar;	// 送信バイナリデータのバイト数
			break;
		case UDP_UINT16_T:
			ret = sizeof(uint16_t)*NumOfVar;// 送信バイナリデータのバイト数
			break;
		case UDP_UINT32_T:
			ret = sizeof(uint32_t)*NumOfVar;// 送信バイナリデータのバイト数
			break;
		case UDP_DOUBLE:
			ret = sizeof(double)*NumOfVar;	// 送信バイナリデータのバイト数
			break;
	}
	
	return ret;
}

void UDPTransmitter::SetNumOfVar(unsigned int NumOfVar){
	// 送信する変数の数を再設定する関数
	// 注意：コンストラクタのNumOfVarを超えてはならない！
	TransBytes = CalcTransBytes(VariableType, NumOfVar);	// 変数の数と送信するデータの型からバイト数を再計算
}

// オーバーロードここから (そのうちテンプレートで書きなおす予定)
void UDPTransmitter::Transmit(const char* Data, unsigned int Length){
	// charデータを送信する関数
	// Data : 送信データ, Length : 送信データのバイト数
	sendto(sock, Data, (size_t)Length, 0,(struct sockaddr *)&addr_dest, sizeof(addr_dest));
}
// 下記のmemcpyはunion使った方が速い？
volatile bool UDPTransmitter::Transmit(int16_t* Data){
	// int16_tデータを送信する関数
	// Data : 送信データ
	memcpy(TransData,(int16_t*)Data,TransBytes);	// charにコピー(エンディアンに注意！)
	return SendTo();	// 送信
}
volatile bool UDPTransmitter::Transmit(int32_t* Data){
	// int32_tデータを送信する関数
	// Data : 送信データ
	memcpy(TransData,(int32_t*)Data,TransBytes);	// charにコピー(エンディアンに注意！)
	return SendTo();	// 送信
}
volatile bool UDPTransmitter::Transmit(uint16_t* Data){
	// uint16_tデータを送信する関数
	// Data : 送信データ
	memcpy(TransData,(uint16_t*)Data,TransBytes);	// charにコピー(エンディアンに注意！)
	return SendTo();	// 送信
}
volatile bool UDPTransmitter::Transmit(uint32_t* Data){
	// uint32_tデータを送信する関数
	// Data : 送信データ
	memcpy(TransData,(uint32_t*)Data,TransBytes);	// charにコピー(エンディアンに注意！)
	return SendTo();	// 送信
}
volatile bool UDPTransmitter::Transmit(double* Data){
	// doubleデータを送信する関数
	// Data : 送信データ
	memcpy(TransData,(double*)Data,TransBytes);		// charにコピー(エンディアンに注意！)
	return SendTo();	// 送信
}
// オーバーロードここまで

