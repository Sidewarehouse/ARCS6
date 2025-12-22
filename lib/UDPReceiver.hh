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

#ifndef UDPRECV
#define UDPRECV

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

namespace ARCS {	// ARCS名前空間
	//! @brief UDP受信器クラス
	class UDPReceiver {
		public:
			// 受信データの型の設定用定義
			enum VarType {
				UDP_INT16_T,	// int16_t型
				UDP_INT32_T,	// int32_t型
				UDP_UINT16_T,	// uint16_t型
				UDP_UINT32_T,	// uint32_t型
				UDP_DOUBLE		// double型
			};
			// コンストラクタ
			// PortNum : 待受Port番号, TypeOfVar : 受信する変数の型(ヘッダの定義を参照せよ)
			// NumOfVar : 受信する変数の数
			UDPReceiver(const unsigned short PortNum, VarType TypeOfVar, unsigned int NumOfVar);
			UDPReceiver(const char* IPaddress, const unsigned short PortNum, VarType TypeOfVar, unsigned int NumOfVar);
			~UDPReceiver();						// デストラクタ
			volatile bool GetReceiveFlag(void);	// 受信フラグを取得する関数
			volatile bool ReBind(void);			// バインドし直す関数
			// 変数値と実際に受信した要素数を取得する関数(オーバーロード)
			unsigned int Receive(int16_t* Data);
			unsigned int Receive(int32_t* Data);
			unsigned int Receive(uint16_t* Data);
			unsigned int Receive(uint32_t* Data);
			unsigned int Receive(double* Data);
		
		private:
			UDPReceiver(const UDPReceiver&);					// コピーコンストラクタ使用禁止
			const UDPReceiver& operator=(const UDPReceiver&);	// 代入演算子使用禁止
			static const unsigned int TimeToReBind = 100000;	// [us] バインド再試行周期
			static void ReceiveThread(UDPReceiver *p);			// 受信スレッド
			int sock;					// ソケットファイルディスクリプタ
			struct sockaddr_in addr;	// ソケット構造体
			VarType VariableType;		// 受信するデータの型
			unsigned int VarNum;		// 受信する変数の数
			char* RecvData;				// 受信バイナリデータ
			size_t RecvBytes;			// 受信バイナリデータのバイト数
			ssize_t ActRecvBytes;		// 実際に受信したバイナリデータのバイト数(若しくはエラー番号)
			pthread_t ThreadID;			// スレッドID
			volatile bool ReceiveFlag;	// 受信フラグ 受信時にtrue，未受信時にfalseとなる。リセットはGetReceiveFlagから行う。
			static size_t CalcRecvBytes(VarType TypeOfVar, unsigned int NumOfVar);	// 受信するデータの型からバイト数を計算する関数
		
	};
}

#endif



