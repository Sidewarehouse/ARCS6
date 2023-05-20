//! @file pcie_ac01.c
//! @brief PCIe-AC01用デバイスドライバ
//! @date 2017/01/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2017 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the Dual BSD/GPL License.

#include <linux/kernel.h> 
#include <linux/module.h>
#include <linux/kthread.h>
#include <linux/sched.h>
#include <linux/delay.h>
#include <linux/proc_fs.h>
#include <linux/cdev.h>
#include <asm/uaccess.h>

#define PCIE_AC01_DEV_NAME	"pcie_ac01"			//!< デバイス名
#define PCIE_AC01_VERSION	"REV.1020170130\0"	//!< バージョン情報
#define PCIE_AC01_MEMSIZE	4096				//!< PCIe-AC01のメモリサイズ
#define PCIE_AC01_MAX_AXIS	16					//!< 最大モータ軸数
#define PCIE_AC01_TIMEOUT	200					//!< [μs程度] ビジーウェイト用タイムアウト

// モジュール情報設定
MODULE_DESCRIPTION("ARCS PCIe-AC01 Device Driver");
MODULE_AUTHOR("Yuki YOKOKURA, Siderwarehouse Industries");
MODULE_LICENSE("Dual BSD/GPL");

// モジュールパラメータ設定
static unsigned long MemAddress;			//!< PCIe物理メモリアドレス
module_param(MemAddress, ulong, S_IRUGO);	// insmod引数の設定

// グローバル変数
static dev_t DeviceID;						//!< デバイス番号
static struct cdev DeviceStruct;			//!< キャラクタデバイス構造体
static u32* MemMapping = NULL;				//!< 物理メモリアクセス用マッピング
static volatile bool BusyWaitFlag = false;	//!< 通信待機ビジーウェイトフラグ

// read/write用バッファ

#define MAX_BUFLEN 64
unsigned char cdev_buf[ MAX_BUFLEN ];
//static int cdev_buflen = 0;

//! @brief readデータの定義
enum ReadDataList {
	Read_Status			= 0x0001,	//!< ステータスレジスタ
	Read_DaxisCurrent	= 0x0002,	//!< d軸電流
	Read_QaxisCurrent	= 0x0004,	//!< q軸電流
	Read_EncoderData	= 0x0008,	//!< エンコーダデータ
	Read_EchoBack		= 0x0010,	//!< エコーバック
};

//! @brief writeデータの定義
enum WriteDataList {
	Write_Command			= 0x0001,	//!< コマンドレジスタ
	Write_DaxisCurrentRef	= 0x0002,	//!< d軸電流指令
	Write_QaxisCurrentRef	= 0x0004,	//!< q軸電流指令
	Write_EchoForward		= 0x0010,	//!< エコーフォワード
};

//! @brief 送信変数アドレスの定義
enum TransVarAddr {
	AxisCommand		= 0x0,	//!< コマンド
	DaxisCurrentRef	= 0x1,	//!< d軸電流指令
	QaxisCurrentRef	= 0x2,	//!< q軸電流指令
	EchoForward		= 0xF,	//!< エコーフォワード
};

//! @brief 受信変数アドレスの定義
enum RecvVarAddr {
	AxisStatus		= 0x0,	//!< ステータス
	DaxisCurrent	= 0x1,	//!< d軸電流
	QaxisCurrent	= 0x2,	//!< q軸電流
	EncorderData	= 0x3,	//!< エンコーダデータ
	EchoBack		= 0xF,	//!< エコーバック
};

//! @brief コマンドリストの定義
enum AxisCommandList {
	Req_Status		 = 0x10000000,	//!< ステータス要求
	Req_DaxisCurrent = 0x10000001,	//!< d軸電流要求
	Req_QaxisCurrent = 0x10000002,	//!< q軸電流要求
	Req_EncoderData	 = 0x10000003,	//!< エンコーダデータ要求
	Req_EchoBack	 = 0x1000000F,	//!< エコーバック要求
};

// プロトタイプ宣言
static int __init PCIeAC01_module_init(void);
static void __exit PCIeAC01_module_creanup(void);
static int PCIeAC01_open(struct inode* inode, struct file *fp);
static int PCIeAC01_close(struct inode* inode, struct file* fp);
static ssize_t PCIeAC01_read(struct file* fp, char* buf, size_t count, loff_t* offset);
static ssize_t PCIeAC01_write(struct file* fp, const char* buf, size_t count, loff_t* offset);
static inline void WriteDataRegister(unsigned int AddrNum, volatile u32 WriteData);
static inline volatile u32 ReadDataRegister(unsigned int AddrNum);
static inline void SetCommunicationRef(u16 AxisSetting, enum TransVarAddr VarAddr);
static inline volatile u32 GetCommunicationSts(void);
static inline void ClearCommunicationSts(void);
static inline void WaitTransRecv(void);
static inline void SetTransDataReg(unsigned int Axis, enum TransVarAddr VarAddr, volatile u32 WriteData);
static inline volatile u32 GetRecvDataReg(unsigned int Axis, enum RecvVarAddr VarAddr);
static inline void GetEncorderData(u32* EncData);
static inline void GetQaxisCurrent(u32* Iq);
static inline void SetQaxisCurrentRef(u32* Iqref);

// モジュールの定義
module_init(PCIeAC01_module_init);
module_exit(PCIeAC01_module_creanup);

//! @brief ファイルオペレーション構造体
static struct file_operations DeviceOperations = {
    .owner   = THIS_MODULE,
    .open    = PCIeAC01_open,
    .release = PCIeAC01_close,
    .read    = PCIeAC01_read,
    .write   = PCIeAC01_write,
};

//! @brief モジュール初期化関数
static int __init PCIeAC01_module_init(void){
	printk("ARCS PCIe-AC01 Driver: Starting [%s]\n", PCIE_AC01_VERSION);
	
	// キャラクタデバイス番号の動的取得
	if(alloc_chrdev_region(&DeviceID, 0, 1, PCIE_AC01_DEV_NAME) == -1){
		printk("ARCS PCIe-AC01 Driver: alloc_chrdev_region failed.\n" );	// 取得できなかったときの処理
		return -1;
	}
	
	// キャラクタデバイスの初期化
	cdev_init(&DeviceStruct, &DeviceOperations);
	DeviceStruct.owner = THIS_MODULE;
	
	// キャラクタデバイスの登録
    if(cdev_add( &DeviceStruct, DeviceID, 1) == -1){
        printk( KERN_WARNING "ARCS PCIe-AC01 Driver: cdev_add failed.\n" );	// 登録できなかったときの処理
        return -1;
    }
	
	// 物理メモリアクセスのためのマッピング
	MemMapping = (u32*)ioremap_nocache(MemAddress, PCIE_AC01_MEMSIZE);
	if(MemMapping == NULL){
		printk("ARCS PCIe-AC01 Driver: ioremap_nocache failed.\n");		// マッピングできなかったときの処理
		return -1;
	}
	
	// 情報表示
	printk("ARCS PCIe-AC01 Driver: Major Num. = %d, Minor Num. = %d\n", MAJOR(DeviceID), MINOR(DeviceID));
	printk("ARCS PCIe-AC01 Driver: MemAddress = 0x%08lX, MemSize = %d\n", MemAddress, PCIE_AC01_MEMSIZE);
	printk("ARCS PCIe-AC01 Driver: MemMapping = 0x%p\n", MemMapping);
	
	return 0;
}

//! @brief モジュール消去関数
static void __exit PCIeAC01_module_creanup(void){
	BusyWaitFlag = false;					// 通信待機ビジーウェイトフラグを無効に設定
	iounmap(MemMapping);					// マッピング解除
    cdev_del(&DeviceStruct);				// キャラクタデバイス削除
    unregister_chrdev_region(DeviceID, 1);	// デバイス番号の返却
	printk("ARCS PCIe-AC01 Driver: Ending\n");
}

//! @brief /dev/pcie_ac01をopenしたときに呼ばれる関数
static int PCIeAC01_open(struct inode* inode, struct file* fp){
	printk("ARCS PCIe-AC01 Driver: Device opened.\n");
	ClearCommunicationSts();	// 通信状態レジスタを初期化
	BusyWaitFlag = true;		// 通信待機ビジーウェイトフラグを有効に設定
	return 0;
}

//! @brief /dev/pcie_ac01をcloseしたときに呼ばれる関数
static int PCIeAC01_close(struct inode* inode, struct file* fp){
	BusyWaitFlag = false;		// 通信待機ビジーウェイトフラグを無効に設定
    printk("ARCS PCIe-AC01 Driver: Device closed.\n");
    return 0;
}

//! @brief /dev/pcie_ac01がユーザ空間からreadされたときに呼ばれる関数
//! @param[in]	buf		ユーザ空間上のバッファへのポインタ
//! @param[in]	count	読み込みデータ選択入力 (readデータの定義を参照のこと)
static ssize_t PCIeAC01_read(struct file* fp, char* buf, size_t count, loff_t* offset){
	static u32 ReadDataBuff[PCIE_AC01_MAX_AXIS] = {0};					// ユーザ空間への送信データバッファ
	static const int ReadDataLength = sizeof(u32)*PCIE_AC01_MAX_AXIS;	// ユーザ空間への送信データバッファのサイズ計算
	// 読み込みデータの判別と処理
	switch((int)count) {
		case (int)Read_Status:
			break;
		case (int)Read_DaxisCurrent:
			break;
		case (int)Read_QaxisCurrent:
			GetQaxisCurrent(ReadDataBuff);	// q軸電流の読み込み
			break;
		case (int)Read_EncoderData:
			GetEncorderData(ReadDataBuff);	// エンコーダデータの読み込み
			break;
		case (int)Read_EchoBack:
			break;
		default:
			memset(ReadDataBuff, 0x00, ReadDataLength);	// どれにも当てはまらない場合はバッファクリア
			break;
	}
	// ユーザ空間へのデータ送信
	if(copy_to_user(buf, ReadDataBuff, ReadDataLength)){
		printk("ARCS PCIe-AC01 Driver: copy_to_user failed.\n");	// 送信できなかったときの処理
		return -EFAULT;
	}
	*offset += ReadDataLength;
    return ReadDataLength;
}

//! @brief /dev/pcie_ac01がユーザ空間からwriteされたときに呼ばれる関数
//! @param[in]	buf		ユーザ空間上のバッファへのポインタ
//! @param[in]	count	書き込みデータ選択入力 (writeデータの定義を参照のこと)
static ssize_t PCIeAC01_write(struct file* fp, const char* buf, size_t count, loff_t* offset){
	static u32 WriteDataBuff[PCIE_AC01_MAX_AXIS] = {0};					// ユーザ空間からの受信データバッファ
	static const int WriteDataLength = sizeof(u32)*PCIE_AC01_MAX_AXIS;	// ユーザ空間からの受信データバッファのサイズ計算
	// ユーザ空間からデータ受信
	if(copy_from_user(WriteDataBuff, buf, WriteDataLength)){
		printk("ARCS PCIe-AC01 Driver: copy_from_user failed.\n");	// 受信できなかったときの処理
		return -EFAULT;
	}
	// 書き込みデータの判別と処理
	switch((int)count) {
		case (int)Write_Command:
			break;
		case (int)Write_DaxisCurrentRef:
			break;
		case (int)Write_QaxisCurrentRef:
			SetQaxisCurrentRef(WriteDataBuff);	// q軸電流指令の書き込み
			break;
		case (int)Write_EchoForward:
			break;
		default:
			break;
	}
	//printk("!!!!!! WriteDataBuff = %d, count = %d\n", WriteDataBuff[0], (int)count);
	*offset += WriteDataLength;
	return WriteDataLength;
}

//! @brief データレジスタに書き込む関数
//! @param[in]	AddrNum		アドレス要素番号(ベースアドレスから32bit区切りで何個目か)
//! @param[in]	WriteData	書き込みデータ
static inline void WriteDataRegister(unsigned int AddrNum, volatile u32 WriteData){
	MemMapping[AddrNum] = WriteData;	// 書き込み実行
}

//! @brief データレジスタから読み込む関数
//! @param[in]	AddrNum		アドレス要素番号(ベースアドレスから32bit区切りで何個目か)
//! @return 読み込みデータ
static inline volatile u32 ReadDataRegister(unsigned int AddrNum){
	return MemMapping[AddrNum];	// 読み込み実行
}

//! @brief 通信指令レジスタへ書き込む関数
//! @param[in]	AxisSetting	軸設定ビットパターン
//! @param[in]	VarAddr		送信変数アドレス
static inline void SetCommunicationRef(u16 AxisSetting, enum TransVarAddr VarAddr){
	volatile u32 WriteData = (((volatile u32)AxisSetting) << 16) + VarAddr;	// 軸設定と送信変数アドレスから書き込むべきデータを計算
	WriteDataRegister(0, WriteData);
}

//! @brief 通信状態レジスタから読み込む関数
//! @return 通信状態レジスタの値
static inline volatile u32 GetCommunicationSts(void){
	return ReadDataRegister(1);
}

//! @brief 通信状態レジスタをクリアする関数
static inline void ClearCommunicationSts(void){
	WriteDataRegister(1, 0x00000000);
}

//! @brief データ送信/受信を待機する関数(ブロッキング関数)
static inline void WaitTransRecv(void){
	// 通信状態レジスタビジーウェイト
	unsigned int TimeOut = PCIE_AC01_TIMEOUT;	// タイムアウト設定
	while(BusyWaitFlag){
		if(GetCommunicationSts() == 0x00000001 || TimeOut == 0) break;	// 送信/受信を完了 OR タイムアウトしたら無限ループを抜ける
		TimeOut--;			// タイムアウト減算
		usleep_range(1,2);	// 1～2μsの間スリープする
	}
	ClearCommunicationSts();	// 通信状態レジスタをクリア
}

//! @brief 送信データレジスタへ書き込む関数
//! @param[in]	Axis		軸番号(1から始まる)
//! @param[in]	VarAddr		送信変数アドレス
//! @param[in]	WriteData	送信変数データ
static inline void SetTransDataReg(unsigned int Axis, enum TransVarAddr VarAddr, volatile u32 WriteData){
	unsigned int OfstAddrNum = (Axis - 1)*32 + 4;		// 軸番号から先頭アドレス番号を計算
	WriteDataRegister(OfstAddrNum + VarAddr, WriteData);// 指定された送信データレジスタに書き込む
}

//! @brief 受信データレジスタから読み込む関数
//! @param[in]	Axis		軸番号(1から始まる)
//! @param[in]	VarAddr		受信変数アドレス
//! @return		受信変数データ
static inline volatile u32 GetRecvDataReg(unsigned int Axis, enum RecvVarAddr VarAddr){
	unsigned int OfstAddrNum = (Axis - 1)*32 + 16 + 4;	// 軸番号から先頭アドレス番号を計算
	return ReadDataRegister(OfstAddrNum + VarAddr);		// 指定された受信データレジスタから読み込む
}

//! @brief エンコーダデータを読み込む関数
static inline void GetEncorderData(u32* EncData){
	SetTransDataReg(1, AxisCommand, Req_EncoderData);	// コマンドにエンコーダデータ要求を設定 1軸目
	SetCommunicationRef(0x0001, AxisCommand);			// コマンド送信指令
	WaitTransRecv();									// 受信完了待機
	EncData[0] = GetRecvDataReg(1, EncorderData);		// エンコーダデータレジスタからデータ取得 1軸目
}

//! @brief q軸電流を読み込む関数
static inline void GetQaxisCurrent(u32* Iq){
	SetTransDataReg(1, AxisCommand, Req_QaxisCurrent);	// コマンドにq軸電流要求を設定 1軸目
	SetCommunicationRef(0x0001, AxisCommand);			// コマンド送信指令
	WaitTransRecv();									// 受信完了待機
	Iq[0] = GetRecvDataReg(1, QaxisCurrent);			// q軸電流レジスタからデータ取得 1軸目
}

//! @brief q軸電流指令を書き込む関数
static inline void SetQaxisCurrentRef(u32* Iqref){
	SetTransDataReg(1, QaxisCurrentRef, Iqref[0]);		// q軸電流指令を設定 1軸目
	SetCommunicationRef(0x0001, QaxisCurrentRef);		// q軸電流指令送信指令
	WaitTransRecv();									// 送信完了待機
}

