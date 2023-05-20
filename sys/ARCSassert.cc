//! @file ARCSassert.cc
//! @brief ARCS用ASSERTクラス
//! @date 2020/03/12
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <unistd.h>
#include <stdio.h>
#include <tuple>
#include "ARCSassert.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

// 静的メンバ変数の実体
pthread_mutex_t ARCSassert::SyncMutex = PTHREAD_MUTEX_INITIALIZER;	//!< Mutexロック
pthread_cond_t ARCSassert::SyncCond = PTHREAD_COND_INITIALIZER;		//!< 待機条件
bool ARCSassert::RealtimeModeFlag = false;		//!< リアルタイムモードフラグ true = リアルタイム，false = 非リアルタイム
bool ARCSassert::EmergencyStopFlag = false;		//!< 緊急停止フラグ  true = 緊急停止, false = 通常
bool ARCSassert::EmergencyProcDoneFlag = true;	//!< 緊急停止処理完了フラグ true = 完了，false = 未完了 (WaitEmergencyが呼ばれる前はtrueで初期化しておく，WaitEmergencyが呼ばれる前にassertが失敗するとassert_from_macroでブロックされてしまうのを防止)
std::string ARCSassert::EmergencyStopCond = "";	//!< 緊急停止時該当条件
std::string ARCSassert::EmergencyStopFile = "";	//!< 緊急停止時ファイル名
int ARCSassert::EmergencyStopLine = 0;			//!< 緊急停止時行番号

//! @brief コンストラクタ
ARCSassert::ARCSassert(void){
	PassedLog();
	pthread_mutex_init(&SyncMutex, nullptr);	// Mutex初期化
	pthread_cond_init(&SyncCond, nullptr);		// 条件初期化
}

//! @brief デストラクタ
ARCSassert::~ARCSassert(){
	PassedLog();
}

//! @brief ARCS用assert関数(マクロから呼ばれることを想定)
//! @param[in] assertion assert条件
//! @param[in] condition assert条件の文字列
//! @param[in] file assert条件に引っかかったときのファイル名
//! @param[in] line assert条件に引っかかったときの行番号
void ARCSassert::assert_from_macro(const bool assertion, const std::string& condition, const std::string& file, const int line){
	if(assertion == false && EmergencyStopFlag == false){	// assertion条件に一致しなくて且つ以前に緊急停止されていないとき
		// assert情報を設定
		pthread_mutex_lock(&SyncMutex);		// Mutexロック
		EmergencyStopFlag = true;			// 緊急停止フラグを立てる
		EmergencyStopCond = condition;		// assert条件の格納
		EmergencyStopFile = file;			// assertion failed のときのファイル名の格納
		EmergencyStopLine = line;			// assertion failed のときの行番号の格納
		pthread_cond_broadcast(&SyncCond);	// assertion情報が更新されたことを配信
		pthread_mutex_unlock(&SyncMutex);	// Mutexロック解除
		PassedLog();
		
		// 緊急停止処理が完了するまで待機
		EventLog("Waiting for Emergency Process...");
		pthread_mutex_lock(&SyncMutex);		// Mutexロック
		while(EmergencyProcDoneFlag != true){
			pthread_cond_wait(&SyncCond, &SyncMutex);	// 緊急停止処理が完了するまで待機
		}
		pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
		EventLog("Waiting for Emergency Process...Done");
		
		// イベントログに状況を書き込む
		ARCSeventlog::WriteEventLog(condition, ">>>>> EMERGENCY STOP <<<<< : " + file, line);
		
		if(RealtimeModeFlag == true){
			// リアルタイムモードで動作中のとき
			EventLog("Infinity Blocking at Line of Assertion Failed......");
			while(1) usleep(WAIT_TIME);		// Assertion Failedのコード位置から先には進ませない
		}else{
			// 非リアルタイムモードで動作中のとき
			// assertionメッセージ表示
			PassedLog();
			printf("ASSERTION FAILED IN NON-REALTIME MODE\n");			// エラーメッセージ表示
			printf("  CONDITION   : %s\n", EmergencyStopCond.c_str());	// 引っ掛かった条件表示
			printf("  FILE NAME   : %s\n", EmergencyStopFile.c_str());	// 引っ掛かったファイル名
			printf("  LINE NUMBER : %d\n", EmergencyStopLine);			// 引っ掛かった行番号
			exit(1);													// 強制終了
		}
	}
}

//! @brief リアルタイムモードかどうか返す関数
//! @return フラグ true = リアルタイム，false = 非リアルタイム
bool ARCSassert::IsRealtimeMode(void){
	bool ret;
	pthread_mutex_lock(&SyncMutex);
	ret = RealtimeModeFlag;
	pthread_mutex_unlock(&SyncMutex);
	return ret;
}

//! @brief リアルタイムモードに切り替える関数
void ARCSassert::SetRealtimeMode(void){
	pthread_mutex_lock(&SyncMutex);
	RealtimeModeFlag = true;
	pthread_mutex_unlock(&SyncMutex);
	EventLog("Switched to Realtime Mode of Emergency Stop.");
}

//! @brief 非リアルタイムモードに切り替える関数
void ARCSassert::SetNonRealtimeMode(void){
	pthread_mutex_lock(&SyncMutex);
	RealtimeModeFlag = false;
	pthread_mutex_unlock(&SyncMutex);
	EventLog("Switched to Non-Realtime Mode of Emergency Stop.");
}

//! @brief 緊急事態かどうか返す関数
//! @return フラグ true = 緊急事態，false = 正常
bool ARCSassert::IsEmergency(void){
	bool ret;
	pthread_mutex_lock(&SyncMutex);
	ret = EmergencyStopFlag;
	pthread_mutex_unlock(&SyncMutex);
	return ret;
}

//! @brief 緊急事態を宣言する関数
void ARCSassert::DeclareEmergency(void){
	pthread_mutex_lock(&SyncMutex);
	EmergencyStopFlag = true;	// PANPANPAN, this is ARCS60, we declare emergency!
	pthread_mutex_unlock(&SyncMutex);
}

//! @brief 緊急処理が完了したことを知らせる関数
void ARCSassert::DoneEmergencyProc(void){
	PassedLog();
	pthread_mutex_lock(&SyncMutex);
	EmergencyProcDoneFlag = true;
	pthread_cond_broadcast(&SyncCond);
	pthread_mutex_unlock(&SyncMutex);
	PassedLog();
}

//! @brief 緊急停止処理が完了したかどうか返す関数
//! @return フラグ true = 完了，false = 未完了
bool ARCSassert::IsDoneEmergencyProc(void){
	bool ret;
	pthread_mutex_lock(&SyncMutex);
	ret = EmergencyProcDoneFlag;
	pthread_mutex_unlock(&SyncMutex);
	return ret;
}

//! @brief Assert条件，ファイル名，行番号を返す関数
//! @return 条件，ファイル名，行番号
std::tuple< std::string, std::string, int > ARCSassert::GetAssertInfo(void){
	std::string cond, file;
	int line;
	pthread_mutex_lock(&SyncMutex);
	cond = EmergencyStopCond;
	file = EmergencyStopFile;
	line = EmergencyStopLine;
	pthread_mutex_unlock(&SyncMutex);
	return std::forward_as_tuple(cond, file, line);
}

//! @brief 緊急停止が発生するまで待機する関数
void ARCSassert::WaitEmergency(void){
	EmergencyProcDoneFlag = false;		// 初期化(WaitEmergencyが呼ばれる前にassertが失敗するとassert_from_macroでブロックされてしまうのを防止)
	EventLog("Waiting for Emergency...");
	pthread_mutex_lock(&SyncMutex);		// Mutexロック
	while(EmergencyStopFlag != true){
		pthread_cond_wait(&SyncCond, &SyncMutex);	// 緊急停止が発生するまで待機
	}
	pthread_mutex_unlock(&SyncMutex);	// Mutexアンロック
	EventLog("Waiting for Emergency...Done");
}


