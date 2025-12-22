//! @file UDPReceiver.hh
//! @brief UDP受信器クラス
//!
//! 指定したIPアドレスおよびPort番号からUDPで配列変数値を受け取るクラス
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
#include <pthread.h>
#include "UDPReceiver.hh"

using namespace ARCS;

UDPReceiver::UDPReceiver(const unsigned short PortNum, VarType TypeOfVar, unsigned int NumOfVar)
	// コンストラクタ (関数オーバーロード：すべてのIPアドレスを指定する場合)
	// PortNum : 待受Port番号, TypeOfVar : 受信する変数の型(ヘッダの定義を参照せよ)
	// NumOfVar : 受信する変数の数
: sock(0), addr(), VariableType(TypeOfVar), VarNum(NumOfVar), RecvData(0), RecvBytes(0),
  ActRecvBytes(0), ThreadID(0), ReceiveFlag(false)
{
	sock = socket(AF_INET, SOCK_DGRAM, 0);	// ソケット生成
	addr.sin_family = AF_INET;				// ファミリ指定
	addr.sin_port = htons(PortNum);			// Port番号の設定
	addr.sin_addr.s_addr = INADDR_ANY;		// すべてのIPアドレス
	
	RecvBytes = CalcRecvBytes(TypeOfVar, NumOfVar);	// 受信するデータの型からバイト数を計算する
	RecvData = new char[RecvBytes];					// 受信バイナリデータ用の配列確保
	memset(RecvData, 0, RecvBytes);					// 配列初期化
	
	struct sched_param scp;	// スレッドパラメータ
	scp.sched_priority = 1;	// スレッド優先度を1に設定
	pthread_create(&ThreadID, NULL, (void*(*)(void*))ReceiveThread, this);	// 受信スレッド生成＆受信開始
	pthread_setschedparam(ThreadID, SCHED_FIFO, &scp);	// スレッド優先度設定
}

UDPReceiver::UDPReceiver(const char* IPaddress, const unsigned short PortNum, VarType TypeOfVar, unsigned int NumOfVar)
	// コンストラクタ (関数オーバーロード：固有のIPアドレスを指定した場合)
	// PortNum : 待受Port番号, TypeOfVar : 受信する変数の型(ヘッダの定義を参照せよ)
	// NumOfVar : 受信する変数の数
: sock(0), addr(), VariableType(TypeOfVar), VarNum(NumOfVar), RecvData(0), RecvBytes(0),
  ActRecvBytes(0), ThreadID(0), ReceiveFlag(false)
{
	sock = socket(AF_INET, SOCK_DGRAM, 0);	// ソケット生成
	addr.sin_family = AF_INET;				// ファミリ指定
	addr.sin_port = htons(PortNum);			// Port番号の設定
	addr.sin_addr.s_addr = inet_addr(IPaddress);	// IPアドレスの設定
	
	RecvBytes = CalcRecvBytes(TypeOfVar, NumOfVar);	// 受信するデータの型からバイト数を計算する
	RecvData = new char[RecvBytes];					// 受信バイナリデータ用の配列確保
	memset(RecvData, 0, RecvBytes);					// 配列初期化
	
	struct sched_param scp;	// スレッドパラメータ
	scp.sched_priority = 1;	// スレッド優先度を1に設定
	pthread_create(&ThreadID, NULL, (void*(*)(void*))ReceiveThread, this);	// 受信スレッド生成＆受信開始
	pthread_setschedparam(ThreadID, SCHED_FIFO, &scp);	// スレッド優先度設定
}


UDPReceiver::~UDPReceiver(){
	// デストラクタ
	// ソケット破棄
	pthread_cancel(ThreadID);		// スレッド終了
	pthread_join(ThreadID, NULL);	// スレッド解放
	close(sock);
	delete[] RecvData;	// 受信バイナリデータ用配列の破棄
}

size_t UDPReceiver::CalcRecvBytes(VarType TypeOfVar, unsigned int NumOfVar){
	// 受信するデータの型からバイト数を計算する関数 (そのうち気が向いたらテンプレートで書きなおすかも)
	// (わざわざこんなことをしているのは、このクラスを使う人がバイト数を意識しないで済むようにするため)
	
	size_t ret = 0;
	
	switch(TypeOfVar){
		case UDP_INT16_T:
			ret = sizeof(int16_t)*NumOfVar;	// 受信バイナリデータのバイト数
			break;
		case UDP_INT32_T:
			ret = sizeof(int32_t)*NumOfVar;	// 受信バイナリデータのバイト数
			break;
		case UDP_UINT16_T:
			ret = sizeof(uint16_t)*NumOfVar;// 受信バイナリデータのバイト数
			break;
		case UDP_UINT32_T:
			ret = sizeof(uint32_t)*NumOfVar;// 受信バイナリデータのバイト数
			break;
		case UDP_DOUBLE:
			ret = sizeof(double)*NumOfVar;	// 受信バイナリデータのバイト数
			break;
	}
	
	return ret;
}

void UDPReceiver::ReceiveThread(UDPReceiver *p){
	// 受信スレッド
	
	// バインドできるまで無限ループ
	while(1){
		// バインドに成功したら無限ループを抜ける
		if(bind(p->sock, (struct sockaddr *)&(p->addr), sizeof(p->addr)) == 0)break;
		usleep(p->TimeToReBind);	// 次のバインド再試行まで待機
	}
	
	// 受信のための無限ループ
	while(1){
		p->ActRecvBytes = recv(p->sock, p->RecvData, p->RecvBytes, 0);	// UDP受信(ブロッキングモード)
		p->ReceiveFlag = true;							// 受信フラグをセット
	}
}

volatile bool UDPReceiver::GetReceiveFlag(void){
	// 受信フラグを取得する関数
	volatile bool ret = ReceiveFlag;
	ReceiveFlag = false;	// 受信フラグのリセット
	return ret;				// 受信フラグを返す
}

volatile bool UDPReceiver::ReBind(void){
	// バインドし直す関数
	// バインドに成功したらtrueを返す
	if(bind(sock, (struct sockaddr *)&addr, sizeof(addr)) == 0)return true;
	// 失敗のときはfalseを返す
	return false;
}

// 関数オーバーロードここから (そのうち気が向いたらテンプレートで書きなおすかも)
unsigned int UDPReceiver::Receive(int16_t* Data){
	// int16_t型の変数値と実際に受信した要素数を取得する関数
	memcpy(Data,RecvData,RecvBytes);
	return ActRecvBytes/sizeof(int16_t);
}
unsigned int UDPReceiver::Receive(int32_t* Data){
	// int32_t型の変数値と実際に受信した要素数を取得する関数
	memcpy(Data,RecvData,RecvBytes);
	return ActRecvBytes/sizeof(int32_t);
}
unsigned int UDPReceiver::Receive(uint16_t* Data){
	// uint16_t型の変数値と実際に受信した要素数を取得する関数
	memcpy(Data,RecvData,RecvBytes);
	return ActRecvBytes/sizeof(uint16_t);
}
unsigned int UDPReceiver::Receive(uint32_t* Data){
	// uint32_t型の変数値と実際に受信した要素数を取得する関数
	memcpy(Data,RecvData,RecvBytes);
	return ActRecvBytes/sizeof(uint32_t);
}
unsigned int UDPReceiver::Receive(double* Data){
	// double型の変数値と実際に受信した要素数を取得する関数
	memcpy(Data,RecvData,RecvBytes);
	return ActRecvBytes/sizeof(double);
}
// 関数オーバーロードここまで

