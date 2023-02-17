//! @file ARCScommon.cc
//! @brief ARCS共通静的関数クラス
//!
//! ARCS内で共通に使用する静的関数クラス
//!
//! @date 2020/04/07
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <unistd.h>
#include <time.h>
#include <pthread.h>
#include "ARCScommon.hh"
#include "ARCSeventlog.hh"
#include "ARCSassert.hh"

using namespace ARCS;

//! @brief スレッドを動作させるCPUとポリシーを設定する関数
//! @param[in]	FuncName	スレッド関数ID
//! @param[in]	CPUnum		CPUコア番号(0スタート)
//! @param[in]	Policy		SCHED_FIFO とか SCHED_RR などのポリシー
//! @param[in]	Priority	優先順位(高い 0 ←→ 99 低い)
void ARCScommon::SetCPUandPolicy(const pthread_t FuncName, const unsigned int CPUnum, const int Policy, const int Priority){
	// CPUコアの割り当て
	cpu_set_t cpuset;			// CPU設定用変数
	CPU_ZERO(&cpuset);			// まずCPU設定用変数をクリアして，
	CPU_SET(CPUnum, &cpuset);	// 指定したCPUコア番号をCPU設定用変数にセットし，
	arcs_assert(pthread_setaffinity_np(FuncName, sizeof(cpu_set_t), &cpuset) == 0);	// スレッドをCPUコアに割り当てる
	arcs_assert(pthread_getaffinity_np(FuncName, sizeof(cpu_set_t), &cpuset) == 0);	// 実際にスレッドが，
	arcs_assert(CPU_ISSET(CPUnum, &cpuset) == 1);	// 指定したCPUコア番号に設定されているかをチェック
	
	// 優先度設定
	struct sched_param FuncName_param;							// スレッドパラメータ
	FuncName_param.sched_priority = sched_get_priority_max(Policy) - Priority;	// 可能な限り高い優先度から優先度を設定(0～99)
	pthread_setschedparam(FuncName, Policy, &FuncName_param);	// 優先度を可能な限り高く設定
	
	// イベントログ
	EventLogVar(CPUnum);
	EventLogVar(Policy);
	EventLogVar(FuncName_param.sched_priority);
}

//! @brief 浮動小数点から文字列へ (std::to_stringだとフォーマット指定できないぽいので仕方なく実装)
//! @param[in] u 浮動小数点数
//! @param[in] format 表示フォーマット (printfと同一の書式)
//! @return 変換後の文字列
std::string ARCScommon::DoubleToString(double u, const std::string& format){
	char charbuff[NUM_CHARBUFF] = {'\0'};	// 数値→文字列変換バッファ
	sprintf(charbuff, format.c_str(), u);	// 文字列へ変換
	return std::string(charbuff);
}

//! @brief 整数値uint64_tから文字列へ
//! @param[in] u uint64_t
//! @param[in] format 表示フォーマット (printfと同一の書式)
//! @return 変換後の文字列
std::string ARCScommon::Uint64ToString(uint64_t u, const std::string& format){
	char charbuff[NUM_CHARBUFF] = {'\0'};	// 数値→文字列変換バッファ
	sprintf(charbuff, format.c_str(), u);	// 文字列へ変換
	return std::string(charbuff);
}

//! @brief 現在時刻を取得する関数
//! @return 現在時刻(文字列)
std::string ARCScommon::GetNowTime(void){
	time_t NowTime;
	time(&NowTime);	// 現在時刻の取得
	return ctime(&NowTime);
}

