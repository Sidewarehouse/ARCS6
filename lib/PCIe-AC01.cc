//! @file PCIe-AC01.cc
//! @brief PCIe-AC01 インタフェースクラス
//! @date 2018/12/30
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2018 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#include <cassert>
#include <fcntl.h>
#include <pthread.h>
#include "PCIe-AC01.hh"

using namespace ARCS;

//! @brief コンストラクタ
//! @param[in]	BaseAddress	PCIe-AC01ボードのベースアドレス
//! @param[in]	CPUNum		デバイスドライバを駆動するCPU番号
PCIeAC01::PCIeAC01(const off_t BaseAddress, unsigned int CPUNum)
	: ShellPath("/bin/sh"),							// シェルまでのファイルパス
	  StartPath("./mod/pcie_ac01/StartModule.sh"),	// 開始スクリプトまでのファイルパス
	  StopPath("./mod/pcie_ac01/StopModule.sh"),	// 停止スクリプトまでのファイルパス
	  fd(0), cpuset()
{
	// ドライバモジュール開始
	char buff[11];
	sprintf(buff, "0x%08X", (uint32_t)BaseAddress);
	std::string AddrParam = buff;
	ExecLinuxCommand(ShellPath + " " + StartPath + " " + AddrParam);
	
	// CPUコア割り当て用
	CPU_ZERO(&cpuset);			// CPU設定用変数をクリアして，
	CPU_SET(CPUNum, &cpuset);	// 所望のCPUコア番号をCPU設定用変数にセット
	
	// ファイルを開く
	fd = open("/dev/pcie_ac01", O_RDWR | O_NOATIME);
	assert(0 <= fd && "/dev/pcie_ac01 cannot be opened.\n");	// 開けなかったときの処理
}

//! @brief デストラクタ
PCIeAC01::~PCIeAC01(){
	close(fd);	// ファイルを閉じる
	ExecLinuxCommand(ShellPath + " " + StopPath);	// ドライバモジュール停止
}

//! @brief Linuxコマンドを実行する関数
//! @param[in]	CommandLine	Linuxコマンド
void PCIeAC01::ExecLinuxCommand(std::string CommandLine){
	FILE* fp = popen(CommandLine.c_str(), "w");
	assert(fp != nullptr);
	pclose(fp);
}

//! @brief エンコーダデータを読み込む関数
//! @param[out]	EncData	エンコーダデータ
void PCIeAC01::GetEncorderData(uint32_t EncData[MaxAxis]){
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);	// スレッドをCPUコアに割り当てる
	ssize_t ret = read(fd, EncData, (size_t)Read_EncoderData);	// エンコーダ受信
	assert(ret != 0);
}

//! @brief q軸電流を読み込む関数
//! @param[out]	Iq	q軸電流
void PCIeAC01::GetQaxisCurrent(double Iq[MaxAxis]){
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);	// スレッドをCPUコアに割り当てる
	uint32_t_float Iq_buff;	// uint32_t→float変換用共用体
	ssize_t ret = read(fd, Iq_buff.binary_data, (size_t)Read_QaxisCurrent);	// q軸電流受信
	assert(ret != 0);
	for(unsigned int i=0;i<MaxAxis;i++){
		Iq[i] = (double)Iq_buff.float_data[i];	// float→doubleのキャスト
	}
}

//! @brief q軸電流指令を書き込む関数
//! @param[in]	Iqref	q軸電流指令
void PCIeAC01::SetQaxisCurrentRef(double Iqref[MaxAxis]){
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);	// スレッドをCPUコアに割り当てる
	uint32_t_float Iq_buff;	// float→uint32_t変換用共用体
	for(unsigned int i=0;i<MaxAxis;i++){
		Iq_buff.float_data[i] = (float)Iqref[i];	// double→floatのキャスト
	}
	ssize_t ret = write(fd, Iq_buff.binary_data, (size_t)Write_QaxisCurrentRef);	// q軸電流指令送信
	assert(ret != -1);
}

