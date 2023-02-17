//! @file AmpIncSquareWave.hh
//! @brief 振幅増加方形波生成器
//!
//! 時間に従って振幅が増加していく方形波を生成する関数
//!
//! @date 2021/12/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include "AmpIncSquareWave.hh"
#include "SquareWave.hh"
#include "StairsWave.hh"

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

//! @brief 振幅増加方形波生成器
//! @param[in]	Time	[s]  時刻
//! @param[in]	Freq	[Hz] 周波数
//! @param[in]	Ystp	[-/period]  方形波1周期毎の振幅増加量
//! @param[in]	Nstp	[-] 方形波周期の数
//! @return	出力
double ARCS::AmpIncSquareWave(const double Time, const double Freq, const double Ystp, const double Nstp){
	double y1 = SquareWave(Freq, 0, Time);					// 方形波生成
	double y2 = StairsWave(Time, 0, Ystp, 1.0/Freq, Nstp);	// 階段波生成
	return y1*y2;
}
