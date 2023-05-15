//! @file CPUSettings.cc
//! @brief CPU設定クラス
//!
//! CPUのコア，ポリシーなどを設定するクラス
//!
//! @date 2019/07/25
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "CPUSettings.hh"
#include <cassert>

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

using namespace ARCS;

//! @brief スレッドを動作させるCPUコアとポリシーを設定する関数
//! @param[in]	FuncName	スレッド関数ID
//! @param[in]	CPUnum		CPUコア番号(0スタート)
//! @param[in]	Policy		SCHED_FIFO とか SCHED_RR などのポリシー
void CPUSettings::SetCPUandPolicy(const pthread_t FuncName, const unsigned int CPUnum, const int Policy){
	// CPUコアの割り当て
	cpu_set_t cpuset;			// CPU設定用変数
	CPU_ZERO(&cpuset);			// まずCPU設定用変数をクリアして，
	CPU_SET(CPUnum, &cpuset);	// 指定したCPUコア番号をCPU設定用変数にセットし，
	arcs_assert(pthread_setaffinity_np(FuncName, sizeof(cpu_set_t), &cpuset) == 0);	// スレッドをCPUコアに割り当てる
	arcs_assert(pthread_getaffinity_np(FuncName, sizeof(cpu_set_t), &cpuset) == 0);	// 実際にスレッドが，
	arcs_assert(CPU_ISSET(CPUnum, &cpuset) == 1);	// 指定したCPUコア番号に設定されているかをチェック
	
	// 優先度設定
	struct sched_param FuncName_param;								// スレッドパラメータ
	FuncName_param.sched_priority = sched_get_priority_max(Policy);	// 可能な限り高い優先度を取得(0～99)
	pthread_setschedparam(FuncName, Policy, &FuncName_param);		// 優先度を可能な限り高く設定
	
	// イベントログ
	EventLogVar(CPUnum);
	EventLogVar(Policy);
	EventLogVar(FuncName_param.sched_priority);
}

