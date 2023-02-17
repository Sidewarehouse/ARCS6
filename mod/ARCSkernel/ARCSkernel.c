//! @file ARCSkernel.c
//! @brief ARCSリアルタイムカーネルモジュール
//! @date 2019/08/01
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
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

#define ARCSKERNEL_DEV_NAME	"arcs_kernel"	//!< デバイス名
#define ARCSKERNEL_VERSION	"REV.19080115"	//!< バージョン情報
#define ARCSKERNEL_MEMSIZE	4096			//!< PCIe-AC01のメモリサイズ
#define ARCSKERNEL_WRDATNUM	16				//!< 書き込みデータの数

// モジュール情報設定
MODULE_DESCRIPTION("ARCS Realtime Kernel Module");
MODULE_AUTHOR("Yuki YOKOKURA, Siderwarehouse Laboratory");
MODULE_LICENSE("Dual BSD/GPL");

// モジュールパラメータ設定
//static unsigned long MemAddress;			//!< PCIe物理メモリアドレス
//module_param(MemAddress, ulong, S_IRUGO);	// insmod引数の設定

// グローバル変数
static dev_t DeviceID;						//!< デバイス番号
static struct cdev DeviceStruct;			//!< キャラクタデバイス構造体
static u32* MemMapping = NULL;				//!< 物理メモリアクセス用マッピング
static volatile bool BusyWaitFlag = false;	//!< 通信待機ビジーウェイトフラグ

// read/write用バッファ
#define MAX_BUFLEN 64
unsigned char cdev_buf[ MAX_BUFLEN ];
//static int cdev_buflen = 0;

// プロトタイプ宣言
static int __init ARCSkernel_module_init(void);
static void __exit ARCSkernel_module_creanup(void);
static int ARCSkernel_open(struct inode* inode, struct file *fp);
static int ARCSkernel_close(struct inode* inode, struct file* fp);
static ssize_t ARCSkernel_read(struct file* fp, char* buf, size_t count, loff_t* offset);
static ssize_t ARCSkernel_write(struct file* fp, const char* buf, size_t count, loff_t* offset);

// モジュールの定義
module_init(ARCSkernel_module_init);
module_exit(ARCSkernel_module_creanup);

//! @brief ファイルオペレーション構造体
static struct file_operations DeviceOperations = {
    .owner   = THIS_MODULE,
    .open    = ARCSkernel_open,
    .release = ARCSkernel_close,
    .read    = ARCSkernel_read,
    .write   = ARCSkernel_write,
};

//! @brief モジュール初期化関数
static int __init ARCSkernel_module_init(void){
	printk("ARCS Realtime Kernel Module: Starting [%s]\n", ARCSKERNEL_VERSION);
	
	// キャラクタデバイス番号の動的取得
	if(alloc_chrdev_region(&DeviceID, 0, 1, ARCSKERNEL_DEV_NAME) == -1){
		printk("ARCS Realtime Kernel Module: alloc_chrdev_region failed.\n" );	// 取得できなかったときの処理
		return -1;
	}
	
	// キャラクタデバイスの初期化
	cdev_init(&DeviceStruct, &DeviceOperations);
	DeviceStruct.owner = THIS_MODULE;
	
	// キャラクタデバイスの登録
    if(cdev_add( &DeviceStruct, DeviceID, 1) == -1){
        printk( KERN_WARNING "ARCS Realtime Kernel Module: cdev_add failed.\n" );	// 登録できなかったときの処理
        return -1;
    }
	
	// 物理メモリアクセスのためのマッピング
	MemMapping = (u32*)ioremap_nocache(MemAddress, ARCSKERNEL_MEMSIZE);
	if(MemMapping == NULL){
		printk("ARCS Realtime Kernel Module: ioremap_nocache failed.\n");		// マッピングできなかったときの処理
		return -1;
	}
	
	// 情報表示
	printk("ARCS Realtime Kernel Module: Major Num. = %d, Minor Num. = %d\n", MAJOR(DeviceID), MINOR(DeviceID));
	printk("ARCS Realtime Kernel Module: MemAddress = 0x%08lX, MemSize = %d\n", MemAddress, ARCSKERNEL_MEMSIZE);
	printk("ARCS Realtime Kernel Module: MemMapping = 0x%p\n", MemMapping);
	
	return 0;
}

//! @brief モジュール消去関数
static void __exit ARCSkernel_module_creanup(void){
	BusyWaitFlag = false;					// 通信待機ビジーウェイトフラグを無効に設定
	iounmap(MemMapping);					// マッピング解除
    cdev_del(&DeviceStruct);				// キャラクタデバイス削除
    unregister_chrdev_region(DeviceID, 1);	// デバイス番号の返却
	printk("ARCS Realtime Kernel Module: Ending\n");
}

//! @brief /dev/arcs_kernelをopenしたときに呼ばれる関数
static int ARCSkernel_open(struct inode* inode, struct file* fp){
	printk("ARCS Realtime Kernel Module: Device opened.\n");
	return 0;
}

//! @brief /dev/arcs_kernelをcloseしたときに呼ばれる関数
static int ARCSkernel_close(struct inode* inode, struct file* fp){
    printk("ARCS Realtime Kernel Module: Device closed.\n");
    return 0;
}

//! @brief /dev/arcs_kernelがユーザ空間からreadされたときに呼ばれる関数
//! @param[in]	buf		ユーザ空間上のバッファへのポインタ
//! @param[in]	count	読み込みデータ選択入力 (readデータの定義を参照のこと)
static ssize_t ARCSkernel_read(struct file* fp, char* buf, size_t count, loff_t* offset){
	static u32 ReadDataBuff[ARCSKERNEL_WRDATNUM] = {0};					// ユーザ空間への送信データバッファ
	static const int ReadDataLength = sizeof(u32)*ARCSKERNEL_WRDATNUM;	// ユーザ空間への送信データバッファのサイズ計算
	// ユーザ空間へのデータ送信
	if(copy_to_user(buf, ReadDataBuff, ReadDataLength)){
		printk("ARCS Realtime Kernel Module: copy_to_user failed.\n");	// 送信できなかったときの処理
		return -EFAULT;
	}
	*offset += ReadDataLength;
    return ReadDataLength;
}

//! @brief /dev/arcs_kernelがユーザ空間からwriteされたときに呼ばれる関数
//! @param[in]	buf		ユーザ空間上のバッファへのポインタ
//! @param[in]	count	書き込みデータ選択入力 (writeデータの定義を参照のこと)
static ssize_t ARCSkernel_write(struct file* fp, const char* buf, size_t count, loff_t* offset){
	static u32 WriteDataBuff[ARCSKERNEL_WRDATNUM] = {0};					// ユーザ空間からの受信データバッファ
	static const int WriteDataLength = sizeof(u32)*ARCSKERNEL_WRDATNUM;	// ユーザ空間からの受信データバッファのサイズ計算
	// ユーザ空間からデータ受信
	if(copy_from_user(WriteDataBuff, buf, WriteDataLength)){
		printk("ARCS Realtime Kernel Module: copy_from_user failed.\n");	// 受信できなかったときの処理
		return -EFAULT;
	}
	*offset += WriteDataLength;
	return WriteDataLength;
}

