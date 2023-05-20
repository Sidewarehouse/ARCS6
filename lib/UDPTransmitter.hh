//! @file UDPTransmitter.hh
//! @brief UDP送信器クラス
//!
//! 指定したIPアドレスおよびPort番号へUDPで配列変数値を飛ばすクラス
//!
//! @date 2019/07/19
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef UDPTRANS
#define UDPTRANS

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

namespace ARCS {	// ARCS名前空間
	//! @brief UDP送信器クラス
	class UDPTransmitter {
		public:
			// 送信データの型の設定用定義
			enum VarType {
				UDP_INT16_T,	// int16_t型
				UDP_INT32_T,	// int32_t型
				UDP_UINT16_T,	// uint16_t型
				UDP_UINT32_T,	// uint32_t型
				UDP_DOUBLE		// double型
			};
			// コンストラクタ (関数オーバーロード: 宛先のみ指定する場合)
			// IPaddress_dest : 宛先IPアドレス, PortNum_dest : 宛先Port番号
			// TypeOfVar : 送信する変数の型(ヘッダの定義を参照せよ), NumOfVar : 送信する変数の数
			UDPTransmitter(const char* IPaddress, const unsigned short PortNum, VarType TypeOfVar, unsigned int NumOfVar);
			// コンストラクタ (関数オーバーロード: 発信元を指定する場合)
			// IPaddress_dept : 発信元IPアドレス, PortNum_dept : 発信元Port番号
			// IPaddress_dest : 宛先IPアドレス, PortNum_dest : 宛先Port番号
			// TypeOfVar : 送信する変数の型, NumOfVar : 送信する変数の数
			UDPTransmitter(const char* IPaddress_dept, const unsigned short PortNum_dept,
				const char* IPaddress_dest, const unsigned short PortNum_dest, VarType TypeOfVar, unsigned int NumOfVar);
			~UDPTransmitter();	// デストラクタ
			void SetNumOfVar(unsigned int NumOfVar);	// 送信する変数の数を再設定する関数  注意：コンストラクタのNumOfVarを超えてはならない！
			// オーバーロードここから (そのうちテンプレートで書きなおす予定)
			// charデータを送信する関数
			// Data : 送信データ, Length : 送信データのバイト数
			void Transmit(const char* Data, unsigned int Length);
			// int16_tデータを送信する関数 Data : 送信データ
			volatile bool Transmit(int16_t* Data);
			// int32_tデータを送信する関数 Data : 送信データ
			volatile bool Transmit(int32_t* Data);
			// uint16_tデータを送信する関数 Data : 送信データ
			volatile bool Transmit(uint16_t* Data);
			// uint32_tデータを送信する関数 Data : 送信データ
			volatile bool Transmit(uint32_t* Data);
			// doubleデータを送信する関数 Data : 送信データ
			volatile bool Transmit(double* Data);
			// オーバーロードここまで
		
		private:
			UDPTransmitter(const UDPTransmitter&);					// コピーコンストラクタ使用禁止
			const UDPTransmitter& operator=(const UDPTransmitter&);	// 代入演算子使用禁止
			volatile bool SendTo(void);	// 実際に送信を実行する関数
			int sock;					// ソケットファイルディスクリプタ
			struct sockaddr_in addr_dest;	// 宛先ソケット構造体
			struct sockaddr_in addr_dept;	// 発信元ソケット構造体
			VarType VariableType;		// 送信するデータの型
			char* TransData;			// 送信バイナリデータ
			size_t TransBytes;			// 送信バイナリデータのバイト数
			static size_t CalcTransBytes(VarType TypeOfVar, unsigned int NumOfVar);	// 送信するデータの型からバイト数を計算する関数
	};
}

#endif

